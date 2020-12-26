#
# Generates sample paths for large community (S=100)
#

using Parameters, 
	Statistics, 
	Random, 
	LinearAlgebra, 
	Distributions, 
	DifferentialEquations, 
	StatsBase, 
	StatsPlots, 
	Plots, 
	DataFrames, 
	CSV, 
	Plotly,
	Alert

pyplot()

include("/home/bb/Gits/white.noise.community.ecology/sde_functions.jl")

# duration of simulation
DURATION = 1e5;

tspan = (0.0,DURATION);

# background parameters
S = 100;
λ = fill(1.0, S);  # niche breadths
U = fill(1.0, S);  # total niche use
c = zeros(S,S);    # strengths of competition
for i in 1:S
	for j in 1:S
		c[i,j] = 1e-6
	end
end
E = fill(0.0, S);  # segregation variances
μ = fill(5.0, S);  # mutation rates
V = fill(1.0, S);  # magnitudes of drift
R = fill(10.0, S);  # innate rate of growth
a = fill(2e-3,S);  # strengths of abiotic selection
θ = fill(0.0, S);  # phenotypic optima

pars = SDEModelParameters(S=S, λ=λ, U=U, E=E, c=c, a=a, μ=μ, V=V, R=R, θ=θ);

#
# find deterministic equilibrium
#

# initial condition
x̄₀ = rand(MvNormal(θ,1),1);
G₀ = rand(LogNormal(4,0.5),S);
N₀ = rand(LogNormal(15,1),S);
u₀ = vec(cat(x̄₀,G₀,N₀,dims=1));

# LogNormal fyi:
#  - parameterization: LogNormal(μ,σ)
#  - expectation: exp( μ + σ²/2 )
#  - variance:    [exp(σ²)-1]exp( 2μ + σ² )

# numerically solve SDE
prob = SDEProblem(fsde,gsde,u₀,tspan,pars);
sol = solve(prob)

TIME = sol.t;
stp = Int64(floor(length(TIME)/500));
TME = 1:stp:length(TIME);

# extract state variables
x̄_sol = transpose(copy(sol[(0*S+1):(1*S),TME]));
G_sol = transpose(copy(sol[(1*S+1):(2*S),TME]));
N_sol = transpose(copy(sol[(2*S+1):(3*S),TME]));

for t in 1:length(TME)
	for i in 1:S
		if N_sol[t,i] < 0.0
			N_sol[t,i] = 0.0
			G_sol[t,i] = 0.0
			x̄_sol[t,i] = Inf
		end
	end
end

gps = fill("Extant",S)
extant_spp  = findall(x->x>0,N_sol[length(TME),:])
extinct_spp = findall(x->x==0,N_sol[length(TME),:])

for spp in extinct_spp
	gps[spp] = "Extinct"
end

# p1 = Plots.plot(TIME[TME],x̄_sol[:,extinct_spp],leg = false,c=:red,linewidth=0.5,ylab="Mean Trait")
# Plots.plot!(TIME[TME],x̄_sol[:,extant_spp],leg = false,c=:black,linewidth=0.5,ylab="Mean Trait")

# p2 = Plots.plot(TIME[TME],log.(G_sol[:,extinct_spp]),leg = false,c=:red,linewidth=0.5,ylab="log(Trait Variance)")
# Plots.plot!(TIME[TME],log.(G_sol[:,extant_spp]),leg = false,c=:black,linewidth=0.5,ylab="log(Trait Variance)")

# p3 = Plots.plot(TIME[TME],log.(N_sol[:,extinct_spp]),leg = false,c=:red,linewidth=0.5,ylab="log(Abundance)",xlab="Time")
# Plots.plot!(TIME[TME],log.(N_sol[:,extant_spp]),leg = false,c=:black,linewidth=0.5,ylab="log(Abundance)",xlab="Time")

# Plots.plot(p1, p2, p3, layout = (3, 1), size=(400,750))

# build dataframes
spp = string("Species", 1)
last = findall(x->x==0,N_sol[:,1])
if length(last) > 0 
	last = last[1] - 1
else
	last = length(TME)
end
df = DataFrame(spp = spp, x = x̄_sol[1:last,1], G = G_sol[1:last,1], N = N_sol[1:last,1], time = TIME[TME[1:last]], ext = gps[1])
for i in 2:S
	spp = string("Species", i)
	last = findall(x->x==0,N_sol[:,i])
	if length(last) > 0 
		last = last[1] - 1
	else
		last = length(TME)
	end
	df_dummy = DataFrame(spp = spp, x = x̄_sol[1:last,i], G = G_sol[1:last,i], N = N_sol[1:last,i], time = TIME[TME[1:last]], ext = gps[i])
	global df = append!(df,df_dummy)
end

# weak competition: c = 1e-8
# CSV.write("/home/bb/Gits/white.noise.community.ecology/sample_path_wc.csv", df)

# moderate competition: c = 1e-6
CSV.write("/home/bb/Gits/white.noise.community.ecology/sample_path_mc.csv", df)

alert("DONE!")