################################################################################
##
## AUTHOR: Bob Week
##
## DATE: 10/16/2020
##
## In this script we provide definitions of data structures for model parameters
## and for state variables used in our simulations of individual-based models.
##
################################################################################

# data type that holds population and parameters

@with_kw mutable struct population
    x::Vector{Float64}  # trait values
    LT::Vector{Float64} # life times (drawn iid from Exp(1) distr)
    N::Int64            # poulation size
    x̄::Float64          # mean trait
    σ²::Float64         # phenotypic variance
    R::Float64          # innate rate of growth
    a::Float64          # strength of abiotic selection
    θ::Float64          # abiotic optimum
    c::Float64          # strength of competition
    μ::Float64          # rate of diffusion (mutation)
    V::Float64          # variance in reproductive output
end

@with_kw mutable struct community
    S::Int64                   # number of species
    x::Vector{Vector{Float64}} # trait values
    LT::Vector{Float64}        # life times (drawn iid from Exp(1) distr)
    g::Vector{Vector{Float64}} # breeding values
    N::Vector{Int64}           # population sizes
    n::Vector{Int64}           # index of rescaling
    x̄::Vector{Float64}         # mean traits
    σ²::Vector{Float64}        # phenotypic variances
    G::Vector{Float64}         # additive genetic variances
    R::Vector{Float64}         # innate rates of growth
    a::Vector{Float64}         # strengths of abiotic selection
    θ::Vector{Float64}         # abiotic optima
    c::Vector{Float64}         # strengths of competition
    λ::Vector{Float64}         # individual niche widths
    U::Vector{Float64}         # total niche uses
    η::Vector{Float64}         # segregation variance
    μ::Vector{Float64}         # rates of diffusion (mutation)
    V::Vector{Float64}         # variances in reproductive output
end

#
# The below methods are for incrementing the discrete-time population processes
# corresponding to non-overlapping generations. Mathematically, such a
# process is considered an extension of a branching random walk.
#
# We provide two versions of the discrete-time update method. Both versions
# extend the basic model of a branching random walk to include competition for
# resources and abiotic stabilizing selection on a quantitative trait.
#
# The first version ("single_indep") tracks just a single species (it ignores
# species i>1) and assumes competition between individuals is independent of
# their trait values.
#
# The second version ("comm_update") tracks S species and assumes competition
# between individuals rapidly diminishes as their trait values diverge.
#

# update for single species
function disc_single_indep(X)

    @unpack S, x, g, N, n, x̄, σ², G, R, a, θ, c, λ, U, η, μ, V = X

    #
    x̄ₚ = fill(0.0,S)
    σₚ²= fill(0.0,S)
    Gₚ = fill(0.0,S)
    Nₚ = fill(0.0,S)

    # creates array of offspring
    # breeding and trait values
    # first index is species
    # second index is individual
    gₚ = fill(zeros(0),S)
    xₚ = fill(zeros(0),S)

    for i in 1:S

        w = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            #

            w̄ = exp( ( R[i] - (a[i]*(θ[i]-x[i][j])^2/2.0) - c[i]*N[i]/n[i] ) / n[i] )

            # parameterizing the NegativeBinomial
            q = w̄/V[i]
            s = w̄^2/(V[i]-w̄)

            # draw random number of offspring
            w[j] = rand( NegativeBinomial( s, q ), 1)[1]

        end

        # tracks the current offspring
        count = Int64(1)

        # loop through parents
        for j in 1:N[i]

            # birth each offspring
            for k in 1:w[j]

                # draw random breeding value for this individual
                append!( gₚ[i], rand( Normal( g[i][j], √(μ[i]/n[i]) ), 1)[1] )

                # draw random trait value for this individual
                append!( xₚ[i], rand( Normal( gₚ[i][count], √η[i] ), 1)[1] )

                count += 1

            end

        end

        x̄ₚ[i] = mean(xₚ[i])
        σₚ²[i]= var(xₚ[i])
        Gₚ[i] = var(gₚ[i])
        Nₚ[i] = sum(w)

    end


    Xₚ = community(S=S,x=xₚ,g=gₚ,N=Nₚ,n=n,x̄=x̄ₚ,σ²=σₚ²,G=Gₚ,R=R,
        a=a,θ=θ,c=c,λ=λ,U=U,η=η,μ=μ,V=V)

    return Xₚ

end

# update for community with non-overlapping generations
function disc_comm_update(X)

    @unpack S, x, N, x̄, σ², R, a, θ, c, λ, U, μ, V = X

    # creates array of offspring trait values
    # first index is species
    # second index is individual
    xₚ = fill(zeros(0),S)

    for i in 1:S

        w = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            #

            # container for aggregating effects of competition
            B = 0.0

            # collect effects of competition with other individuals
            # within the same population
            for k in filter(x -> x≠j, 1:N[i])
                B += U[i]^2*exp( (x[i][j] - x[i][k])^2 / (4*λ[i]) ) / √(4*π*λ[i])
            end

            # collect effects of competition with other individuals
            # in other populations
            for k in filter(x -> x≠i, 1:S)
                for l in 1:N[k]
                    B += U[i]*U[k]*exp( (x[i][j] - x[k][l])^2 / (2*(λ[i]+λ[k])) ) / √(2*π*(λ[i]+λ[k]))
                end
            end

            w̄ = exp( R[i] - a[i]*(θ[i]-x[i][j])^2/2.0 - c[i]*B )

            # parameterizing the NegativeBinomial
            q = w̄/V
            s = w̄^2/(V-w̄)

            # draw random number of offspring
            w[j] = rand( NegativeBinomial( s, q ), 1)[1]

        end

        # total number of offspring
        Nₚ = sum(w)

        # container for locations of offspring
        xₚ = fill(0.0,Nₚ)

        # keeps track of which individual is being born
        ct = 0

        # loop throug parents
        for j in 1:N[i]

            # birth each offspring
            for k in 1:W[j]

                # consider next individual
                ct += 1

                # draw random trait for this individual
                xₚ[ct] = rand( Normal( x[i,j], √μ[i] ), 1)[1]

            end

        end

        x̄ₚ[i] = mean(xₚ)
        σₚ²[i]= var(xₚ)

    end

    Xₚ = community(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,R=R,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end

# rescaled update for community with non-overlapping generations
function rescaled_update(X)

    @unpack S, x, N, x̄, σ², R, a, θ, c, λ, U, μ, V = X

    # creates array of offspring trait values
    # first index is species
    # second index is individual
    xₚ = fill(zeros(0),S)

    for i in 1:S

        w = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            #

            # container for aggregating effects of competition
            B = 0.0

            # collect effects of competition with other individuals
            # within the same population
            for k in filter(x -> x≠j, 1:N[i])
                B += U[i]^2*exp( (x[i,j] - x[i,k])^2 / (4*λ[i]) ) / √(4*π*λ[i])
            end

            # collect effects of competition with other individuals
            # in other populations
            for k in filter(x -> x≠i, 1:S)
                for l in 1:N[k]
                    B += U[i]*U[k]*exp( (x[i,j] - x[k,l])^2 / (2*(λ[i]+λ[k])) ) / √(2*π*(λ[i]+λ[k]))
                end
            end

            w̄ = exp( R[i] - a[i]*(θ[i]-x[i,j])^2/2.0 - c[i]*B )

            # parameterizing the NegativeBinomial
            q = w̄/V
            s = w̄^2/(V-w̄)

            # draw random number of offspring
            w[j] = rand( NegativeBinomial( s, q ), 1)[1]

        end

        # total number of offspring
        Nₚ = sum(w)

        # container for locations of offspring
        xₚ = fill(0.0,Nₚ)

        # keeps track of which individual is being born
        ct = 0

        # loop throug parents
        for j in 1:N[i]

            # birth each offspring
            for k in 1:W[j]

                # consider next individual
                ct += 1

                # draw random trait for this individual
                xₚ[ct] = rand( Normal( x[i,j], √μ[i] ), 1)[1]

            end

        end

        x̄ₚ[i] = mean(xₚ)
        σₚ²[i]= var(xₚ)

    end

    Xₚ = community(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,R=R,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end

# update for community with non-overlapping generations
# using lower bound on fitness
function update_lower(X)

    @unpack S, x, N, x̄, σ², R, a, θ, c, w, U, μ, V = X

    # creates array of offspring trait values
    # first index is species
    # second index is individual
    xₚ = fill(zeros(0),S)

    for i in 1:S

        W = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            #

            # container for aggregating effects of competition
            B = 0.0

            # collect effects of competition with other individuals
            # within the same population
            for k in filter(x -> x≠j, 1:N[i])
                B += U[i]^2*exp( (x[i,j] - x[i,k])^2 / (4*w[i]) ) / √(4*π*w[i])
            end

            # collect effects of competition with other individuals
            # in other populations
            for k in filter(x -> x≠i, 1:S)
                for l in 1:N[k]
                    B += U[i]*U[k]*exp( (x[i,j] - x[k,l])^2 / (2*(w[i]+w[k])) ) / √(2*π*(w[i]+w[k]))
                end
            end

            w = exp( R[i] - a[i]*(θ[i]-x[i,j])^2/2.0 - c[i]*B )

            # parameterizing the NegativeBinomial
            q = w/V
            s = w^2/(V-w)

            # draw random number of offspring
            W[j] = rand( NegativeBinomial( s, q ), 1)[1]

        end

        # total number of offspring
        Nₚ = sum(W)

        # container for locations of offspring
        xₚ = fill(0.0,Nₚ)

        # keeps track of which individual is being born
        ct = 0

        # loop throug parents
        for j in 1:N[i]

            # birth each offspring
            for k in 1:W[j]

                # consider next individual
                ct += 1

                # draw random trait for this individual
                xₚ[ct] = rand( Normal( x[i,j], √μ[i] ), 1)[1]

            end

        end

        x̄ₚ[i] = mean(xₚ)
        σₚ²[i]= var(xₚ)

    end

    Xₚ = community(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,R=R,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end


#
# The below methods are for incrementing the continous-time population processes
# corresponding to overlapping generations. Mathematically, such processes can
# be considered extensions of branching Brownian motions.
#
# We provide two versions of the continous-time update method in analogy to
# above.
#

# update for single species
function cont_single_indep(X)

    @unpack S, x, g, N, n, x̄, σ², G, R, a, θ, c, λ, U, η, μ, V = X

    #
    x̄ₚ = fill(0.0,S)
    σₚ²= fill(0.0,S)
    Gₚ = fill(0.0,S)
    Nₚ = fill(0.0,S)

    # creates array of offspring
    # breeding and trait values
    # first index is species
    # second index is individual
    gₚ = fill(zeros(0),S)
    xₚ = fill(zeros(0),S)

    for i in 1:S

        w = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            #

            w̄ = exp( ( R[i] - (a[i]*(θ[i]-x[i][j])^2/2.0) - c[i]*N[i]/n[i] ) / n[i] )

            # parameterizing the NegativeBinomial
            q = w̄/V[i]
            s = w̄^2/(V[i]-w̄)

            # draw random number of offspring
            w[j] = rand( NegativeBinomial( s, q ), 1)[1]

        end

        # tracks the current offspring
        count = Int64(1)

        # loop through parents
        for j in 1:N[i]

            # birth each offspring
            for k in 1:w[j]

                # draw random breeding value for this individual
                append!( gₚ[i], rand( Normal( g[i][j], √(μ[i]/n[i]) ), 1)[1] )

                # draw random trait value for this individual
                append!( xₚ[i], rand( Normal( gₚ[i][count], √η[i] ), 1)[1] )

                count += 1

            end

        end

        x̄ₚ[i] = mean(xₚ[i])
        σₚ²[i]= var(xₚ[i])
        Gₚ[i] = var(gₚ[i])
        Nₚ[i] = sum(w)

    end


    Xₚ = community(S=S,x=xₚ,g=gₚ,N=Nₚ,n=n,x̄=x̄ₚ,σ²=σₚ²,G=Gₚ,R=R,
        a=a,θ=θ,c=c,λ=λ,U=U,η=η,μ=μ,V=V)

    return Xₚ

end

# update for community with non-overlapping generations
function cont_comm_update(X)

    @unpack S, x, N, x̄, σ², R, a, θ, c, λ, U, μ, V = X

    # creates array of offspring trait values
    # first index is species
    # second index is individual
    xₚ = fill(zeros(0),S)

    for i in 1:S

        w = fill(0,N[i])

        for j in 1:N[i]

            #
            # mean fitness of individual j in species i
            #

            # container for aggregating effects of competition
            B = 0.0

            # collect effects of competition with other individuals
            # within the same population
            for k in filter(x -> x≠j, 1:N[i])
                B += U[i]^2*exp( (x[i][j] - x[i][k])^2 / (4*λ[i]) ) / √(4*π*λ[i])
            end

            # collect effects of competition with other individuals
            # in other populations
            for k in filter(x -> x≠i, 1:S)
                for l in 1:N[k]
                    B += U[i]*U[k]*exp( (x[i][j] - x[k][l])^2 / (2*(λ[i]+λ[k])) ) / √(2*π*(λ[i]+λ[k]))
                end
            end

            w̄ = exp( R[i] - a[i]*(θ[i]-x[i][j])^2/2.0 - c[i]*B )

            # parameterizing the NegativeBinomial
            q = w̄/V
            s = w̄^2/(V-w̄)

            # draw random number of offspring
            w[j] = rand( NegativeBinomial( s, q ), 1)[1]

        end

        # total number of offspring
        Nₚ = sum(w)

        # container for locations of offspring
        xₚ = fill(0.0,Nₚ)

        # keeps track of which individual is being born
        ct = 0

        # loop throug parents
        for j in 1:N[i]

            # birth each offspring
            for k in 1:W[j]

                # consider next individual
                ct += 1

                # draw random trait for this individual
                xₚ[ct] = rand( Normal( x[i,j], √μ[i] ), 1)[1]

            end

        end

        x̄ₚ[i] = mean(xₚ)
        σₚ²[i]= var(xₚ)

    end

    Xₚ = community(x=xₚ,N=Nₚ,x̄=x̄ₚ,σ²=σₚ²,R=R,a=a,θ=θ,c=c,μ=μ,V=V)

    return Xₚ

end
