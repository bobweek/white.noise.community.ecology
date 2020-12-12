########################################################################################
##
## AUTHOR: Bob Week
##
## DATE: 10/30/2020
##
## In this script we check our assumption of Gaussian trait distributions for a set of
## competing species. In particular, we solve a system of SPDE where each SPDE 
## corresponds to the abundance density of a species. This script builds directly off
## of 'spde.jl' for solving SPDE.
##
########################################################################################

# load in some libraries
using StochasticDiffEq, 
    LinearAlgebra, 
    Plots, 
    Parameters, 
    PyCall,
    Alert

pyplot()
logocolors = Colors.JULIA_LOGO_COLORS

# defines the resolution of the mesh used to approximate the solution
const X = 200

# defines the mesh
const Mx = Tridiagonal([1.0 for x in 1:X-1],[-2.0 for i in 1:X],[1.0 for i in 1:X-1])

# data type that holds model parameters
@with_kw mutable struct ModelParameters
    μ::Float64	        # mutation rate
    R::Float64      	# intrinsic rates of growth
    a::Float64	        # strength of abiotic selection
    θ::Float64	        # phenotypic optima
    c::Matrix{Float64}  # strength of competition
    λ::Float64          # niche width
    U::Float64          # niche use
    V::Float64	        # variance of reproductive output
    S::Int64            # number of competing species
end

#######################################################
### Define the model
#######################################################

function O(x₁,x₂,U₁,U₂,λ₁,λ₂)
    val = U₁.*U₂.*exp.(-(x₁ .- x₂).^2 ./ (2*(λ₁+λ₂))) ./ √(2*π*(λ₁+λ₂))
    return(val)
end

function f(u,p,t)
    @unpack μ, R, a, θ, c, λ, U, S = p
    du = zeros(S,X)
    for i in 1:S
        du[i,:] = 0.5 .* μ .* Mx*u[i,:] .+ u[i,:].*(R .- 0.5 .* a .* ((1:X) .- θ).^2)
        comp = zeros(X)
        for x in 1:X
            for j in 1:S
                comp[x] += c[i,j]*dot(u[j,:], O(x,1:X,U,U,λ,λ))
            end
        end
        du[i,:] -= u[i,:] .* comp
        for x in 1:X
            if u[i,x]<0           # if the solution becomes negative
                du[i,x] = 0.01    # then quickly restore it towards zero
            end
        end
    end
    return(du)
end

function g(u,p,t)
    @unpack V, S = p
    du = zeros(S,X)
    for i in 1:S
        for x in 1:X
            if u[i,x]>=0
                du[i,x] = √(V*u[i,x])
            end
        end
    end
    return(du)
end

##############################
##############################
##
## SOLVING THE SPDE \/ \/ \/
##
##############################
##############################

###################################################################
### Symmetric Competition, Strong Mutation, Large Population Sizes
###################################################################

# Define the parameter values used
S = 3
μ = 5.0
R = 3.0
a = 0.002
θ = 100 # Niche Locations are shifted by 100 units
c = zeros(S,S)
for i in 1:S
    for j in 1:S
        c[i,j] = 0.001
    end
    c[i,i] = 0.001
end
λ = 1.0
U = 1.0
V = 1.0

DURATION = 50.0

pars = ModelParameters(μ=μ, R=R, a=a, θ=θ, c=c, λ=λ, U=U, V=V, S=S)

# Define the initial condition
ν₀ = zeros(S,X)
x̄₀ = [θ-20,θ,θ+20]
σ₀ = [4, 4, 4]
max₀ = [250.0,250.0,250.0]
for i in 1:S
    supp = (x̄₀[i]-Int64(√σ₀[i])):(x̄₀[i]+Int64(√σ₀[i]))
    for x in supp
        ν₀[i,x] = max₀[i]*exp(1-1.0/(1.0-(x-x̄₀[i])^2/σ₀[i]))
    end
end

prob = SDEProblem(f,g,ν₀,(0.0,DURATION),pars,progress=true,progress_steps=1)
sol = solve(prob,SOSRI())

TIME = sol.t

abun = zeros(S,X,length(TIME))
for t in 1:length(TIME)
    for i in 1:S
        for x in 1:X
            abun[i,x,t] = sol.u[t][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
        end
    end
end

stp = Int64(floor(length(TIME)/10));
TME = 1:stp:length(TIME);


surf1 = plot(TIME[TME],-59:60,abun[1,40:159,TME],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,TME],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,TME],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,TME],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,TME],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,TME],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-sym-comp-stg-mut-lrg-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TIME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TIME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TIME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-sym-comp-stg-mut-lrg-pop.png")

alert("Best Case: Done!")

###################################################################
### Asymmetric Competition, Strong Mutation, Large Population Size
###################################################################

# Define the parameter values used
S = 3
μ = 5.0
R = 3.0
a = 0.002
θ = 100 # Niche Locations are shifted by 100 units
c = zeros(S,S)
for i in 1:S
    for j in 1:S
        c[i,j] = 0.001
    end
    c[i,i] = 0.0005
end
λ = 1.0
U = 1.0
V = 1.0

DURATION = 50.0

pars = ModelParameters(μ=μ, R=R, a=a, θ=θ, c=c, λ=λ, U=U, V=V, S=S)

# Define the initial condition
ν₀ = zeros(S,X)
x̄₀ = [θ-20,θ,θ+20]
σ₀ = [4, 4, 4]
max₀ = [250.0,250.0,250.0]
for i in 1:S
    supp = (x̄₀[i]-Int64(√σ₀[i])):(x̄₀[i]+Int64(√σ₀[i]))
    for x in supp
        ν₀[i,x] = max₀[i]*exp(1-1.0/(1.0-(x-x̄₀[i])^2/σ₀[i]))
    end
end

prob = SDEProblem(f,g,ν₀,(0.0,DURATION),pars)
sol = solve(prob,SOSRI())

TIME = sol.t

abun = zeros(S,X,length(TIME))
for t in 1:length(TIME)
    for i in 1:S
        for x in 1:X
            abun[i,x,t] = sol.u[t][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
        end
    end
end

stp = Int64(floor(length(TIME)/20));
TME = 1:stp:length(TIME);

surf1 = plot(TIME[TME],-59:60,abun[1,40:159,TME],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,TME],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,TME],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,TME],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,TME],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,TME],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-asym-comp-stg-mut-lrg-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TIME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TIME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TIME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-asym-comp-stg-mut-lrg-pop.png")

alert("Asymm Comp: Done!")

################################################################
### Symmetric Competition, Weak Mutation, Large Population Size
################################################################

# Define the parameter values used
S = 3
μ = 0.5
R = 3.0
a = 0.002
θ = 100 # Niche Locations are shifted by 100 units
c = zeros(S,S)
for i in 1:S
    for j in 1:S
        c[i,j] = 0.001
    end
    c[i,i] = 0.001
end
λ = 1.0
U = 1.0
V = 1.0

DURATION = 50.0

pars = ModelParameters(μ=μ, R=R, a=a, θ=θ, c=c, λ=λ, U=U, V=V, S=S)

# Define the initial condition
ν₀ = zeros(S,X)
x̄₀ = [θ-20,θ,θ+20]
σ₀ = [4, 4, 4]
max₀ = [250.0,250.0,250.0]
for i in 1:S
    supp = (x̄₀[i]-Int64(√σ₀[i])):(x̄₀[i]+Int64(√σ₀[i]))
    for x in supp
        ν₀[i,x] = max₀[i]*exp(1-1.0/(1.0-(x-x̄₀[i])^2/σ₀[i]))
    end
end

prob = SDEProblem(f,g,ν₀,(0.0,DURATION),pars)
sol = solve(prob,SOSRI())

TIME = sol.t

abun = zeros(S,X,length(TIME))
for t in 1:length(TIME)
    for i in 1:S
        for x in 1:X
            abun[i,x,t] = sol.u[t][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
        end
    end
end

stp = Int64(floor(length(TIME)/20));
TME = 1:stp:length(TIME);

surf1 = plot(TIME[TME],-59:60,abun[1,40:159,TME],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,TME],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,TME],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,TME],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,TME],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,TME],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-sym-comp-wke-mut-lrg-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TIME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TIME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TIME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-sym-comp-wke-mut-lrg-pop.png")

alert("Weak Mutation: Done!")

####################################################################
### Symmetric Competition, Strong Mutation, Small Population Sizes
####################################################################

# Define the parameter values used
S = 3
μ = 5.0
R = 1.0
a = 0.002
θ = 100 # Niche Locations are shifted by 100 units
c = zeros(S,S)
for i in 1:S
    for j in 1:S
        c[i,j] = 0.05
    end
end
λ = 1.0
U = 1.0
V = 1.0

DURATION = 50.0

pars = ModelParameters(μ=μ, R=R, a=a, θ=θ, c=c, λ=λ, U=U, V=V, S=S)

# Define the initial condition
ν₀ = zeros(S,X)
x̄₀ = [θ-20,θ,θ+20]
σ₀ = [4, 4, 4]
max₀ = [10.0,10.0,10.0]
for i in 1:S
    supp = (x̄₀[i]-Int64(√σ₀[i])):(x̄₀[i]+Int64(√σ₀[i]))
    for x in supp
        ν₀[i,x] = max₀[i]*exp(1-1.0/(1.0-(x-x̄₀[i])^2/σ₀[i]))
    end
end

prob = SDEProblem(f,g,ν₀,(0.0,DURATION),pars)
sol = solve(prob,SOSRI())

TIME = sol.t

abun = zeros(S,X,length(TIME))
for t in 1:length(TIME)
    for i in 1:S
        for x in 1:X
            abun[i,x,t] = sol.u[t][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
        end
    end
end

stp = Int64(floor(length(TIME)/20));
TME = 1:stp:length(TIME);


surf1 = plot(TIME[TME],-59:60,abun[1,40:159,TME],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,TME],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,TME],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,TME],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,TME],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,TME],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-sym-comp-stg-mut-sml-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TIME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TIME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TIME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-sym-comp-stg-mut-sml-pop.png")

alert("Small Pops: Done!")

#
# TODO: 
#   - Redo ind-based sim in cont time (set dt = 0.01 and check when t passes Exp-dist lifetimes)
#   - Solve SPDE system for S = 10 species under Symm Comp, Stg Mut, and Lrg Pops
#   - Track abund, mean & var even when normality fails
#       - Use deterministic PDE and ODE
#