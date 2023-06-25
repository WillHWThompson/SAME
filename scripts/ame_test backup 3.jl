using DrWatson
using Revise
@quickactivate :SAME
using SAME 
using OrdinaryDiffEq, Plots, RecursiveArrayTools
using OrdinaryDiffEq:ODESolution
using GraphPlot
using Plots
N_individuals = 75 
N_groups = 4  
p_dist = Binomial(N_individuals,0.1)                                              
g_dist = Binomial(N_groups,0.5)


p_dist.p
N = 10



# my_binomial(n,u;ϵ = 0.5) = pdf(p_dist,n)*binomial(BigInt(n),BigInt(u))*ϵ^u * (1-ϵ)^(n-u)

# function initialize_u0(N_individuals::Int,p_dist)
#     """
#     initialize_u0 - create an initial distribution for the AME
#     inputs: 
#         N_individuals: the total number of individuals in your population
#         p_dist: the probability distribution used to specify the clique size, from Distributions.jl
#     """
#     G = zeros(N_individuals,N_individuals)
#     G_Inds = CartesianIndices(G)
#     G = map(inds -> G[inds[1],inds[2]] = my_binomial(inds[1],inds[2]),G_Inds)#generate the correct distribution for each occupation number
#     G_norm = G./sum(G)
#     return G_norm
# end


# function voter_model_ame!(dG,G,t,p)
#     N,U  = size(G)
#     for n_prime=1:N,u_prime=1:U
#         u  = u_prime-1
#         n = n_prime-1

#         dG[n_prime,u_prime] = 0
#         dG[n_prime,u_prime] -= G[n_prime,u_prime]*(u-n)*(u/n)#upflip outflux
#         dG[n_prime,u_prime] -= G[n_prime,u_prime]*(u)*((n-u)/n)#downflip outflux
#         u_prime < U && (dG[n_prime,u_prime]+= G[n_prime,u_prime+1]*((u+1)*((n-u-1)/n)))#inwards downflip flux
#         u_prime > 2 && (dG[n_prime,u_prime]+= G[n_prime,u_prime-1]*((n-u+1)*((u-1)/n)))#inwards upflip flux
#     end
#     return dG 
# end


t_max = 500
u₀ = SAME.initialize_u0(N_individuals,p_dist)
dG = zeros(N_individuals,N_individuals)

SAME.voter_model_ame!(dG,u₀,0,0)









p = 0
t_span = (0.,t_max)
prob = ODEProblem(SAME.voter_model_ame!,u₀,t_span,p)
sol = solve(prob,Tsit5(),saveat=1)

vals = map(x -> sum(x),sol)


heatmap(sol[6])
plot!(1:N_individuals,1:N_individuals)
