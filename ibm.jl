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
	StatsBase, StatsPlots, Plots, PlotThemes, DataFrames, CSV, Optim, Alert

pyplot()
theme(:mute)
include("/home/bb/Gits/white.noise.community.ecology/ibm_functions_structs.jl")

########################################################
#
#  An individual-based model for the entire community
#
#  Used to illustrate the diffusion approximation
#
########################################################

# parameter values
λ = 1.0  	# niche breadths
U = 1.0  	# total niche use
c = 1e-5 	# strengths of competition
E = 1e-5 	# environmental variances
μ = 5	 	# mutation rates
V = 1.0  	# magnitudes of drift
R = 0.5  	# innate rate of growth
a = 1e-2 	# strengths of abiotic selection
θ = 0.0  	# phenotypic optima
k = 1.0  	# scaling parameter
n₀ = Int64(100*k) 	# initial number of discrete individuals
N₀ = k*(R-0.5*√(a*μ))/c # initial population mass (set to the equilibrium mass)

# initial breeding values
g₀ = rand(Normal(θ,√√(μ/a)),n₀);

# initial trait values
Eₘ = √E*Matrix(I, n₀, n₀);
x₀ = vec(rand(MvNormal(g₀,Eₘ),1));

# initial lifetimes
LT₀ = [rand(Exponential(1/k),n₀)]

# set up initial population
X = community(S=1, x=fill(x₀,1), g=fill(g₀,1), n₀=fill(n₀,1), k=fill(k,1), N₀=fill(N₀,1), n=fill(n₀,1), x̄=mean.(fill(x₀,1)), σ²=var.(fill(x₀,1)), G=var.(fill(g₀,1)),
	R=fill(R,1), a=fill(a,1), θ=fill(θ,1), c=fill(c,1), λ=fill(λ,1),U=fill(U,1), E=fill(E,1), μ=fill(μ,1), V=fill(V,1), BT=[zeros(n₀)], LT=LT₀ )

CLK = minimum(LT₀[1]) # keeps track of simulation time	

DURATION = 5 # length of time to simulate for

# set up history of population
Xₕ = [X];

# simulate
i = 1
# ppsz = Xₕ[i].n[1]
while CLK < DURATION

	if prod( log.( 1 .+ Xₕ[i].n ) ) > 0

		append!(Xₕ, [cont_single_indep(Xₕ[i])] )

		CLK = minimum( Xₕ[i].LT[1] )

		# ppsz = Xₕ[i].n[1]

		i+=1

	else

		CLK = DURATION

		# ppsz = 100

	end	

end

# traits of each individual across entire simulation
ttl = n₀
inds = zeros(4,ttl)
for i in 1:n₀
	inds[1,i] = Xₕ[1].BT[1][i]
	inds[2,i] = Xₕ[1].LT[1][i]
	inds[3,i] = Xₕ[1].x[1][i]
	inds[4,i] = Xₕ[1].x[1][i]
end
for i in 2:length(Xₕ)
	mnLT = minimum(Xₕ[i-1].LT[1])
	prnt = argmin(Xₕ[i-1].LT[1])
	newbs = findall(x->x==mnLT,Xₕ[i].BT[1])
	ttl += length(newbs)
	indsₚ = zeros(4,length(newbs))
	for j in 1:length(newbs)
		indsₚ[1,j] = Xₕ[i].BT[1][newbs[j]]
		indsₚ[2,j] = Xₕ[i].LT[1][newbs[j]]
		indsₚ[3,j] = Xₕ[i].x[1][newbs[j]]
		indsₚ[4,j] = Xₕ[i-1].x[1][prnt] # parental trait value
	end
	inds = [inds indsₚ]
end

endT = minimum(Xₕ[length(Xₕ)-1].LT[1])
br_pl1 = plot([inds[1,1],inds[2,1]], [inds[3,1],inds[3,1]], c=:black, la=0.7/√k, lw=0.7/√k, legend=false, xrange=[0,endT], ylabel="Trait Value (k=1)", xlabel=("Time"))
for i in 2:ttl
	plot!([inds[1,i],inds[2,i]], [inds[3,i],inds[3,i]], c=:black, la=0.7/√k, lw=0.7/k)
	if inds[3,i]!=inds[4,i]
		plot!([inds[1,i],inds[1,i]], [inds[3,i],inds[4,i]], c=:black, ls=:dot, la=0.7/√k, lw=0.7/√k)
	end
end
br_pl1

rescaled_pl = plot(br_pl10, br_pl5, br_pl1, layout = (3, 1), size=(800,800))

savefig("/home/bb/Gits/white.noise.community.ecology/rescaled_pl.png")

alert("DONE!")

# set up containers for paths of N, x̄ and σ²
nₕ = zeros(1,length(Xₕ))
x̄ₕ = zeros(1,length(Xₕ))
σ²ₕ= zeros(1,length(Xₕ))
Gₕ = zeros(1,length(Xₕ))

# container for individuals
#individualₕ = zeros(2)

# fill them in
for i in 1:1
	for j in 1:length(Xₕ)
		nₕ[i,j] =Xₕ[j].n[i]
		x̄ₕ[i,j] =Xₕ[j].x̄[i]
		σ²ₕ[i,j]=Xₕ[j].σ²[i]
		Gₕ[i,j] =Xₕ[j].G[i]
	end
end


# rescaled time
resc_time = (1:length(Xₕ))./k

plot(resc_time,N₀.*nₕ[1,:]./k)
plot(resc_time,x̄ₕ[1,:]./√k)
plot(resc_time,σ²ₕ[1,:]./k)
