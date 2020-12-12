################################################################################
##
## AUTHOR: Bob Week
##
## DATE: 10/16/2020
##
## In this script we simulate an individual-based model for a single species
## and illustrate what happens in the diffusion limit.
##
## This script depends on another script called "ibm_functions_structs.jl".
## In that script we provide definitions of data structures for model parameters
## and for state variables. In that script we also define methods for iterating
## the simulation. With all the gory details located elsewhere, we can focus
## this script on simulating the model for a specified duration.
##
################################################################################

using Parameters, Statistics, Random, LinearAlgebra, Distributions,
	StatsBase, StatsPlots, Plots, DataFrames, CSV, Optim

include("/home/bb/Gits/white.noise.community.ecology/ibm_functions_structs.jl")

########################################################
#
#  An individual-based model for the entire community
#
#  Used to illustrate the diffusion approximation
#
########################################################

# parameter values
w = 0.1  # niche breadths
U = 1.0  # total niche use
c = 2e-3 # strengths of competition
η = 1e-5 # environmental variances
μ = 1e-3 # mutation rates
V = 2.0  # magnitudes of drift
r = 1.0  # innate rate of growth
a = 1e-2 # strengths of abiotic selection
θ = 0.0  # phenotypic optima
n = 1.0  # scaling parameter

##
## VERY IMPORTANT REQUIREMENT   -->  V >= exp(r)
##
## this inequality must be satisfied to use negative binomial sampling
##

# initial abundance
N₀ = Int64(floor( n*( r -0.5*( η*a + √(μ*a) ) )/c ) )
# here we set this to the equilibrium abundance

# initial breeding values
g₀ = rand(Normal(0.0,1.0),N₀)

# initial trait values
ηₘ = √η*Matrix(I, N₀, N₀)
x₀ = vec(rand(MvNormal(g₀,ηₘ),1))

# set up initial population
X = community(S=1, x=fill(x₀,1), g=fill(g₀,1), N=fill(N₀,1), n=fill(n,1),
	x̄=mean.(fill(x₀,1)), σ²=var.(fill(x₀,1)), G=var.(fill(g₀,1)),
	R=fill(r,1), a=fill(a,1), θ=fill(θ,1), c=fill(c,1), w=fill(w,1),
	U=fill(U,1), η=fill(η,1), μ=fill(μ,1), V=fill(V,1) )

# always a good idea to inspect a single iteration
rescaled_lower(X)

# number of generations to halt at
T = 50

# set up history of population
Xₕ = fill(X,T)

# simulate
for i in 2:T

	if prod( log.( 1 .+ Xₕ[i-1].N ) ) > 0

		Xₕ[i] = rescaled_lower(Xₕ[i-1])

	else

		Xₕ[i] = Xₕ[i-1]

	end

end

# set up containers for paths of N, x̄ and σ²
Nₕ = zeros(1,T)
x̄ₕ = zeros(1,T)
σ²ₕ= zeros(1,T)
Gₕ = zeros(1,T)

# container for individuals
#individualₕ = zeros(2)

# fill them in
for i in 1:1
	for j in 1:T
		Nₕ[i,j] =Xₕ[j].N[i]
		x̄ₕ[i,j] =Xₕ[j].x̄[i]
		σ²ₕ[i,j]=Xₕ[j].σ²[i]
		Gₕ[i,j] =Xₕ[j].G[i]
	end
end

# rescaled time
resc_time = (1:T)./n

# total number of individuals across entire simulation
total_inds = Int64(sum(Nₕ[1,:]))

# traits of each individual across entire simulation
inds = zeros(2,total_inds)

ind = 0
for i in 1:T
	for j in 1:Int64(Nₕ[1,i])

		global ind += 1
		inds[1,ind] = resc_time[i]
		inds[2,ind] = Xₕ[i].x[1][j]

	end
end

scatter(inds[1,:], inds[2,:], legend=false, ms=.5, c=:black)

# build dataframe
df = DataFrame(x = inds[2,:], time = inds[1,:])

CSV.write("/home/bob/Research/Branching Brownian Motion/n_3.csv", df)

plot(resc_time,Nₕ[1,:]./n)
plot(resc_time,x̄ₕ[1,:]./√n)
plot(resc_time,σ²ₕ[1,:]./n)

