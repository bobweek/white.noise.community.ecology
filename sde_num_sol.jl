#
# Generates sample paths for large community (S=100)
#

using Parameters, Statistics, Random, LinearAlgebra, Distributions, DifferentialEquations, StatsBase, StatsPlots, Plots, DataFrames, CSV, Plotly

plotly()

include("/home/bb/Gits/white.noise.community.ecology/sde_functions.jl")

# duration of simulation
DURATION = 1e4;

tspan = (0.0,TIME);

# background parameters
S = 50;
w = fill(0.1, S);  # niche breadths
U = fill(1.0, S);  # total niche use
c = zeros(S,S);      # strengths of competition
for i in 1:S
	for j in 1:S
		c[i,j] = 1e-4
	end
end
Ω = sum(U); # niche use scaling
η = fill(1.0, S);  # segregation variances
μ = fill(1e-1, S); # mutation rates
V = fill(1.0, S);  # magnitudes of drift
R = fill(0.1, S); # innate rate of growth
a = fill(1e-3,S);       # strengths of abiotic selection
θ = fill(0.0, S);  # phenotypic optima

pars = ModelParameters(S=S, w=w, U=U, η=η, c=c, a=a, μ=μ, V=V, R=R, θ=θ, Ω=Ω);

#
# find deterministic equilibrium
#

# initial condition
x̄₀ = rand(MvNormal(θ,1),1);
G₀ = rand(LogNormal(3.6,0.5),S);
N₀ = rand(LogNormal(15,1),S);
u₀ = vec(cat(x̄₀,G₀,N₀,dims=1));

# LogNormal fyi:
#  - parameterization: LogNormal(μ,σ)
#  - expectation: exp( μ + σ²/2 )
#  - variance:    [exp(σ²)-1]exp( 2μ + σ² )

# numerically solve SDE
prob = SDEProblem(f,g,u₀,tspan,pars);
sol = solve(prob)

# extract state variables
x̄_sol = transpose(copy(sol[(0*S+1):(1*S),:]));
G_sol = transpose(copy(sol[(1*S+1):(2*S),:]));
N_sol = transpose(copy(sol[(2*S+1):(3*S),:]));


LT = length(sol);
p1 = Plots.plot(1:LT,x̄_sol,leg = false,ylab="trait value");
p2 = Plots.plot(1:LT,log.(G_sol),leg = false,ylab="log(heritable variation)");
p3 = Plots.plot(1:LT,log.(N_sol),leg = false,ylab="log(abundance)",xlab="time");
Plots.plot(p1, p2, p3, layout = (3, 1), size=(800,1500))
Plots.savefig("/home/bb/Gits/white.noise.community.ecology/pretty")

# build dataframes
spp = string("Species", 1)
df = DataFrame(spp = spp, x = x̄_sol[1,:], G = G_sol[1,:], N = N_sol[1,:], time = sol.t)
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i,:], G = G_sol[i,:],
	N = N_sol[i,:], time = sol.t)
	global df = append!(df,df_dummy)
end

# weak competition
#CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_wc.csv", df)

# moderate competition
#CSV.write("/home/bob/Research/Branching Brownian Motion/sample_path_mc.csv", df)

# for S = 1000

# extract state variables
x̄_sol = copy(.-sol[(0*S+1):(1*S),LT])
G_sol = copy(sol[(1*S+1):(2*S),LT])
N_sol = copy(sol[(2*S+1):(3*S),LT])

# build dataframes
spp = string("Species", 1)
df = DataFrame(spp = spp, x = x̄_sol[1], G = G_sol[1], N = N_sol[1])
for i in 2:S
	spp = string("Species", i)
	df_dummy = DataFrame(spp = spp, x = x̄_sol[i], G = G_sol[i],
	N = N_sol[i])
	global df = append!(df,df_dummy)
end


# weak competition
CSV.write("/home/bb/Gits/branching.brownian.motion.and.spde/sample_path_wc1e3.csv", df)

# moderate competition
#CSV.write("/home/bb/Gits/branching.brownian.motion.and.spde/sample_path_mc1e3.csv", df)
