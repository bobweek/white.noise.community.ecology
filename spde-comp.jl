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
    DifferentialEquations, 
    LinearAlgebra, 
    Plots, 
    Parameters, 
    PyCall,
    Alert

pyplot()
logocolors = Colors.JULIA_LOGO_COLORS

include("/home/bb/Gits/white.noise.community.ecology/sde_functions.jl")

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
        du[i,:] .-= u[i,:] .* comp
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

DURATION = 100.0

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

#
# solve the SPDE
#
prob = SDEProblem(f,g,ν₀,(0.0,DURATION),pars)
sol = solve(prob,SOSRI())

TIME = sol.t

# take a subset of time points
stp = Int64(floor(length(TIME)/30));
TME = 1:stp:length(TIME);

# trait values
Xs = (-X/2+1):(X/2);

# pull out abundance densities, pop sizes, mean traits and trait vars
abun = zeros(S,X,length(TME))
N    = zeros(S,length(TME))
x̄    = zeros(S,length(TME))
σ²   = zeros(S,length(TME))
for t in 1:length(TME)
    for i in 1:S
        for x in 1:X
            abun[i,x,t] = sol.u[TME[t]][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
            N[i,t] = sum(abun[i,:,t])
            x̄[i,t] = dot( Xs, abun[i,:,t] ) / N[i,t]
            σ²[i,t] = dot( (Xs.-x̄[i,t]).^2, abun[i,:,t] ) / N[i,t]
        end
    end
end

surf1 = plot(TIME[TME],-59:60,abun[1,40:159,:],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,:],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,:],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,:],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,:],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,:],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-sym-comp-stg-mut-lrg-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-sym-comp-stg-mut-lrg-pop.png")

#
# solve ODE and compare
#

# initial condition
u₀ = vec(cat(x̄[:,1],σ²[:,1],N[:,1],dims=1))

SDEpars = SDEModelParameters(μ=fill(μ,S), R=fill(R,S), a=fill(a,S), θ=fill(θ-100,S), c=c, λ=fill(λ,S), U=fill(U,S), V=fill(V,S), S=S, E=zeros(S))
ODEprob = ODEProblem(fsde,u₀,(0.0,DURATION),SDEpars);
ODEsol = solve(ODEprob)
x̄_sol = transpose(copy(ODEsol[(0*S+1):(1*S),:]));
G_sol = transpose(copy(ODEsol[(1*S+1):(2*S),:]));
N_sol = transpose(copy(ODEsol[(2*S+1):(3*S),:]));

abunp = plot(TIME[TME],log.(N[1,:]),c=logocolors.red,label="Non-Gaussian",ylabel="Log(Abundance)",xlabel="Time")
plot!(TIME[TME],log.(N[2,:]),c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],log.(N[3,:]),c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,log.(N_sol[:,1]),c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,2]),c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,3]),c=logocolors.green,ls=:dash,label="Gaussian")

mtraitp = plot(TIME[TME],x̄[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Mean Niche\nLocation",xlabel="Time")
plot!(TIME[TME],x̄[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],x̄[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,x̄_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

traitvp = plot(TIME[TME],σ²[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Log(Variance of\nNiche Location)",xlabel="Time")
plot!(TIME[TME],σ²[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],σ²[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,G_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

comparep = plot(abunp, mtraitp, traitvp, layout = (3, 1), size=(800,800))

savefig("/home/bb/Gits/white.noise.community.ecology/compare-sym-comp-stg-mut-lrg-pop.png")

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

DURATION = 100.0

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

# take a subset of time points
stp = Int64(floor(length(TIME)/30));
TME = 1:stp:length(TIME);

# trait values
Xs = (-X/2+1):(X/2);

# pull out abundance densities, pop sizes, mean traits and trait vars
abun = zeros(S,X,length(TME))
N    = zeros(S,length(TME))
x̄    = zeros(S,length(TME))
σ²   = zeros(S,length(TME))
for t in 1:length(TME)
    print("t=",t,"\n")
    for i in 1:S
        print("i=",i,"\n")
        for x in 1:X
            print("x=",x,"\n")
            abun[i,x,t] = sol.u[TME[t]][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
            N[i,t] = sum(abun[i,:,t])
            x̄[i,t] = dot( Xs, abun[i,:,t] ) / N[i,t]
            σ²[i,t] = dot( (Xs.-x̄[i,t]).^2, abun[i,:,t] ) / N[i,t]
        end
    end
end

surf1 = plot(TIME[TME],-59:60,abun[1,40:159,:],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,:],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,:],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,:],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,:],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,:],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-asym-comp-stg-mut-lrg-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-asym-comp-stg-mut-lrg-pop.png")

#
# solve ODE and compare
#

# initial condition
u₀ = vec(cat(x̄[:,1],σ²[:,1],N[:,1],dims=1))

SDEpars = SDEModelParameters(μ=fill(μ,S), R=fill(R,S), a=fill(a,S), θ=fill(θ-100,S), c=c, λ=fill(λ,S), U=fill(U,S), V=fill(V,S), S=S, E=zeros(S))
ODEprob = ODEProblem(fsde,u₀,(0.0,DURATION),SDEpars);
ODEsol = solve(ODEprob)
x̄_sol = transpose(copy(ODEsol[(0*S+1):(1*S),:]));
G_sol = transpose(copy(ODEsol[(1*S+1):(2*S),:]));
N_sol = transpose(copy(ODEsol[(2*S+1):(3*S),:]));

abunp = plot(TIME[TME],log.(N[1,:]),c=logocolors.red,label="Non-Gaussian",ylabel="Log(Abundance)",xlabel="Time")
plot!(TIME[TME],log.(N[2,:]),c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],log.(N[3,:]),c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,log.(N_sol[:,1]),c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,2]),c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,3]),c=logocolors.green,ls=:dash,label="Gaussian")

mtraitp = plot(TIME[TME],x̄[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Mean Niche\nLocation",xlabel="Time")
plot!(TIME[TME],x̄[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],x̄[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,x̄_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

traitvp = plot(TIME[TME],σ²[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Log(Variance of\nNiche Location)",xlabel="Time")
plot!(TIME[TME],σ²[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],σ²[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,G_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

comparep = plot(abunp, mtraitp, traitvp, layout = (3, 1), size=(800,800))

savefig("/home/bb/Gits/white.noise.community.ecology/compare-asym-comp-stg-mut-lrg-pop.png")

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

DURATION = 100.0

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

# take a subset of time points
stp = Int64(floor(length(TIME)/30));
TME = 1:stp:length(TIME);

# trait values
Xs = (-X/2+1):(X/2);

# pull out abundance densities, pop sizes, mean traits and trait vars
abun = zeros(S,X,length(TME))
N    = zeros(S,length(TME))
x̄    = zeros(S,length(TME))
σ²   = zeros(S,length(TME))
for t in 1:length(TME)
    print("t=",t,"\n")
    for i in 1:S
        print("i=",i,"\n")
        for x in 1:X
            print("x=",x,"\n")
            abun[i,x,t] = sol.u[TME[t]][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
            N[i,t] = sum(abun[i,:,t])
            x̄[i,t] = dot( Xs, abun[i,:,t] ) / N[i,t]
            σ²[i,t] = dot( (Xs.-x̄[i,t]).^2, abun[i,:,t] ) / N[i,t]
        end
    end
end

surf1 = plot(TIME[TME],-59:60,abun[1,40:159,:],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,:],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,:],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,:],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,:],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,:],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-sym-comp-wke-mut-lrg-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-sym-comp-wke-mut-lrg-pop.png")

#
# solve ODE and compare
#

# initial condition
u₀ = vec(cat(x̄[:,1],σ²[:,1],N[:,1],dims=1))

SDEpars = SDEModelParameters(μ=fill(μ,S), R=fill(R,S), a=fill(a,S), θ=fill(θ-100,S), c=c, λ=fill(λ,S), U=fill(U,S), V=fill(V,S), S=S, E=zeros(S))
ODEprob = ODEProblem(fsde,u₀,(0.0,DURATION),SDEpars);
ODEsol = solve(ODEprob)
x̄_sol = transpose(copy(ODEsol[(0*S+1):(1*S),:]));
G_sol = transpose(copy(ODEsol[(1*S+1):(2*S),:]));
N_sol = transpose(copy(ODEsol[(2*S+1):(3*S),:]));

abunp = plot(TIME[TME],log.(N[1,:]),c=logocolors.red,label="Non-Gaussian",ylabel="Log(Abundance)",xlabel="Time")
plot!(TIME[TME],log.(N[2,:]),c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],log.(N[3,:]),c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,log.(N_sol[:,1]),c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,2]),c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,3]),c=logocolors.green,ls=:dash,label="Gaussian")

mtraitp = plot(TIME[TME],x̄[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Mean Niche\nLocation",xlabel="Time")
plot!(TIME[TME],x̄[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],x̄[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,x̄_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

traitvp = plot(TIME[TME],σ²[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Log(Variance of\nNiche Location)",xlabel="Time")
plot!(TIME[TME],σ²[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],σ²[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,G_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

comparep = plot(abunp, mtraitp, traitvp, layout = (3, 1), size=(800,800))

savefig("/home/bb/Gits/white.noise.community.ecology/compare-sym-comp-wke-mut-lrg-pop.png")

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

DURATION = 100.0

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

# take a subset of time points
stp = Int64(floor(length(TIME)/30));
TME = 1:stp:length(TIME);

# trait values
Xs = (-X/2+1):(X/2);

# pull out abundance densities, pop sizes, mean traits and trait vars
abun = zeros(S,X,length(TME))
N    = zeros(S,length(TME))
x̄    = zeros(S,length(TME))
σ²   = zeros(S,length(TME))
for t in 1:length(TME)
    print("t=",t,"\n")
    for i in 1:S
        print("i=",i,"\n")
        for x in 1:X
            print("x=",x,"\n")
            abun[i,x,t] = sol.u[TME[t]][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
            N[i,t] = sum(abun[i,:,t])
            x̄[i,t] = dot( Xs, abun[i,:,t] ) / N[i,t]
            σ²[i,t] = dot( (Xs.-x̄[i,t]).^2, abun[i,:,t] ) / N[i,t]
        end
    end
end

surf1 = plot(TIME[TME],-59:60,abun[1,40:159,:],st=:surface,c=cgrad([:black,logocolors.red,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf2 = plot(TIME[TME],-59:60,abun[2,40:159,:],st=:surface,c=cgrad([:black,logocolors.blue,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
surf3 = plot(TIME[TME],-59:60,abun[3,40:159,:],st=:surface,c=cgrad([:black,logocolors.green,:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="Niche Location",xlab="Time",zlab="Abundance\nDensity",zguidefontrotation=90);
cont1 = contour(TIME[TME],-59:60,abun[1,40:159,:],c=cgrad([:black,logocolors.red,:white]),fill=true,ylab="Niche Location",legend=false);
cont2 = contour(TIME[TME],-59:60,abun[2,40:159,:],c=cgrad([:black,logocolors.blue,:white]),fill=true,ylab="Niche Location",legend=false);
cont3 = contour(TIME[TME],-59:60,abun[3,40:159,:],c=cgrad([:black,logocolors.green,:white]),fill=true,ylab="Niche Location",xlab="Time",legend=false);
plot(surf1, cont1, surf2, cont2, surf3, cont3, layout = (3, 2), size=(800,800))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-sym-comp-stg-mut-sml-pop.png")
crv = plot(-59:60,abun[1,40:159,length(TME)],c=logocolors.red,fill=(0,0.2,logocolors.red),legend=false,xlab="Niche Location",ylab="Abundance Density")
plot!(-59:60,abun[2,40:159,length(TME)],c=logocolors.blue,fill=(0,0.2,logocolors.blue),legend=false)
plot!(-59:60,abun[3,40:159,length(TME)],c=logocolors.green,fill=(0,0.2,logocolors.green),legend=false)
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-sym-comp-stg-mut-sml-pop.png")

#
# solve ODE and compare
#

# initial condition
u₀ = vec(cat(x̄[:,1],σ²[:,1],N[:,1],dims=1))

SDEpars = SDEModelParameters(μ=fill(μ,S), R=fill(R,S), a=fill(a,S), θ=fill(θ-100,S), c=c, λ=fill(λ,S), U=fill(U,S), V=fill(V,S), S=S, E=zeros(S))
ODEprob = ODEProblem(fsde,u₀,(0.0,DURATION),SDEpars);
ODEsol = solve(ODEprob)
x̄_sol = transpose(copy(ODEsol[(0*S+1):(1*S),:]));
G_sol = transpose(copy(ODEsol[(1*S+1):(2*S),:]));
N_sol = transpose(copy(ODEsol[(2*S+1):(3*S),:]));

abunp = plot(TIME[TME],log.(N[1,:]),c=logocolors.red,label="Non-Gaussian",ylabel="Log(Abundance)",xlabel="Time")
plot!(TIME[TME],log.(N[2,:]),c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],log.(N[3,:]),c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,log.(N_sol[:,1]),c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,2]),c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,log.(N_sol[:,3]),c=logocolors.green,ls=:dash,label="Gaussian")

mtraitp = plot(TIME[TME],x̄[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Mean Niche\nLocation",xlabel="Time")
plot!(TIME[TME],x̄[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],x̄[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,x̄_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,x̄_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

traitvp = plot(TIME[TME],σ²[1,:],c=logocolors.red,label="Non-Gaussian",ylabel="Log(Variance of\nNiche Location)",xlabel="Time")
plot!(TIME[TME],σ²[2,:],c=logocolors.blue,label="Non-Gaussian")
plot!(TIME[TME],σ²[3,:],c=logocolors.green,label="Non-Gaussian")
plot!(ODEsol.t,G_sol[:,1],c=logocolors.red,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,2],c=logocolors.blue,ls=:dash,label="Gaussian")
plot!(ODEsol.t,G_sol[:,3],c=logocolors.green,ls=:dash,label="Gaussian")

comparep = plot(abunp, mtraitp, traitvp, layout = (3, 1), size=(800,800))

savefig("/home/bb/Gits/white.noise.community.ecology/compare-sym-comp-stg-mut-sml-pop.png")

alert("Small Pops: Done!")

#################################################################################
### Symmetric Competition, Strong Mutation, Large Population Sizes, Ten Species
#################################################################################

# Define the parameter values used
S = 10
μ = 5.0
R = 3.0
a = 0.002
θ = 100 # Niche Locations are shifted by 100 units
c = zeros(S,S)
for i in 1:S
    for j in 1:S
        c[i,j] = 0.0001
    end
end
λ = 1.0
U = 1.0
V = 1.0

DURATION = 100.0

pars = ModelParameters(μ=μ, R=R, a=a, θ=θ, c=c, λ=λ, U=U, V=V, S=S)

# Define the initial condition
ν₀ = zeros(S,X);
x̄₀ = θ.+(-45:10:45);
σ₀ = fill(1,S);
max₀ = fill(250.0,S);
for i in 1:S
    supp = (x̄₀[i]-Int64(√σ₀[i])):(x̄₀[i]+Int64(√σ₀[i]));
    for x in supp
        ν₀[i,x] = max₀[i]*exp(1-1.0/(1.0-(x-x̄₀[i])^2/σ₀[i]));
    end
end

prob = SDEProblem(f,g,ν₀,(0.0,DURATION),pars)
sol = solve(prob,SOSRI())

TIME = sol.t

# take a subset of time points
stp = Int64(floor(length(TIME)/30));
TME = 1:stp:length(TIME);

# trait values
Xs = (-X/2+1):(X/2);

# pull out abundance densities, pop sizes, mean traits and trait vars
abun = zeros(S,X,length(TME))
N    = zeros(S,length(TME))
x̄    = zeros(S,length(TME))
σ²   = zeros(S,length(TME))
for t in 1:length(TME)
    print("t=",t,"\n")
    for i in 1:S
        print("i=",i,"\n")
        for x in 1:X
            print("x=",x,"\n")
            abun[i,x,t] = sol.u[TME[t]][i,x]
            if abun[i,x,t] < 0
                abun[i,x,t] = 0
            end
            N[i,t] = sum(abun[i,:,t])
            x̄[i,t] = dot( Xs, abun[i,:,t] ) / N[i,t]
            σ²[i,t] = dot( (Xs.-x̄[i,t]).^2, abun[i,:,t] ) / N[i,t]
        end
    end
end

clrs = palette(:tab10);
surf = Any[ ]
for i in 1:S
    push!(surf,plot(TIME[TME],-59:60,abun[i,40:159,:],st=:surface,c=cgrad([:black,clrs[i],:white]),linewidth=0.5,linecolor=:black,alpha=0.9,camera=(120,15),legend=false,ylab="",xlab="",zlab="",zguidefontrotation=90))
    push!(surf,contour(TIME[TME],-59:60,abun[i,40:159,:],c=cgrad([:black,clrs[i],:white]),fill=true,ylab="",legend=false))
end
plot(surf..., layout = (5, 4), size=(1200,1100))
savefig("/home/bb/Gits/white.noise.community.ecology/spde-10-spp.png")
crv = plot(-59:60,abun[1,40:159,length(TME)],c=clrs[1],legend=false,xlab="Niche Location",ylab="Abundance Density")
for i in 2:S
    plot!(-59:60,abun[i,40:159,length(TME)],c=clrs[i],legend=false)
end
crv
savefig("/home/bb/Gits/white.noise.community.ecology/eq-dist-10-spp.png")

#
# solve ODE and compare
#

# initial condition
u₀ = vec(cat(x̄[:,1],σ²[:,1],N[:,1],dims=1))

SDEpars = SDEModelParameters(μ=fill(μ,S), R=fill(R,S), a=fill(a,S), θ=fill(θ-100,S), c=c, λ=fill(λ,S), U=fill(U,S), V=fill(V,S), S=S, E=zeros(S))
ODEprob = ODEProblem(fsde,u₀,(0.0,DURATION),SDEpars);
ODEsol = solve(ODEprob)
x̄_sol = transpose(copy(ODEsol[(0*S+1):(1*S),:]));
G_sol = transpose(copy(ODEsol[(1*S+1):(2*S),:]));
N_sol = transpose(copy(ODEsol[(2*S+1):(3*S),:]));

abunp = plot(TIME[TME],log.(N[1,:]),c=clrs[1],label="Non-Gaussian",ylabel="Log(Abundance)",xlabel="Time",legend = false)
for i in 2:S
    plot!(TIME[TME],log.(N[i,:]),c=clrs[i],label="Non-Gaussian")
end
for i in 1:S
    plot!(ODEsol.t,log.(N_sol[:,i]),c=clrs[i],ls=:dash,label="Gaussian")
end
mtraitp = plot(TIME[TME],x̄[1,:],c=clrs[1],label="Non-Gaussian",ylabel="Mean Niche\nLocation",xlabel="Time",legend = false)
for i in 2:S
    plot!(TIME[TME],x̄[i,:],c=clrs[i],label="Non-Gaussian")
end
for i in 1:S
    plot!(ODEsol.t,x̄_sol[:,i],c=clrs[i],ls=:dash,label="Gaussian")
end
traitvp = plot(TIME[TME],σ²[1,:],c=clrs[1],label="Non-Gaussian",ylabel="Log(Variance of\nNiche Location)",xlabel="Time",legend = false)
for i in 2:S
    plot!(TIME[TME],σ²[2,:],c=clrs[i],label="Non-Gaussian")
end
for i in 1:S
    plot!(ODEsol.t,G_sol[:,1],c=clrs[i],ls=:dash,label="Gaussian")
end

comparep = plot(abunp, mtraitp, traitvp, layout = (3, 1), size=(800,800))

savefig("/home/bb/Gits/white.noise.community.ecology/compare-10-spp.png")

alert("Ten Species: Done!")


#
# TODO: 
#   - Redo ind-based sim in cont time (set dt = 0.01 and check when t passes Exp-dist lifetimes)
#   - Solve SPDE system for S = 10 species under Symm Comp, Stg Mut, and Lrg Pops
#   - Track abund, mean & var even when normality fails
#       - Use deterministic PDE and ODE
#