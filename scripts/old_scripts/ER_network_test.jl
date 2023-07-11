
using DrWatson
using Revise
@quickactivate :SAME
using SAME 
#using OrdinaryDiffEq, Plots, RecursiveArrayTools
#using OrdinaryDiffEq:ODESolution
using GraphPlot
using Plots
using Combinatorics
# N_individuals = 75 
# N_groups = 4  
# p_dist = Binomial(N_individuals,0.1)                                              
# g_dist = Binomial(N_groups,0.5)


#N = 10




#p_dist = Binomial(N_individuals,0.1)                                              
#g_dist = Binomial(N_groups,0.5)




# function make_sbm(n_groups,n_individuals,p_dist,g_dist,occupation_prob)
#     edge_list,vertex_attribues = same.generate_hypergraph_bipartite_edge_list(n_groups,n_individuals,p_dist,g_dist)
#     graph,vertex_attribues = same.make_hypergraph(n_groups,n_individuals,p_dist,g_dist) 
#     er_graph = erdos_renyi(n_individuals,occupation_prob)
#     er_edges = edges(er_graph) |> collect
#     graph_combined = union(graph,er_graph)
#     return graph_combined
# end
begin 
N_individuals = 5000 
N_groups = 1000

p_dist = Dirac(N_individuals)                                              
g_dist = Dirac(0)
# edge_list,vertex_attribues = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)
 graph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist) 
graph,hypergraph_meta = SAME.make_SBM(N_groups,N_individuals,p_dist,g_dist,0.01)
#graph2,hypergraph_meta = SAME.make_unipartite_community_graph(N_groups,N_individuals,p_dist,g_dist)
end

# graph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist) 
# edges(graph)
# SAME.bipartite_projection(graph,vertex_attribues,2)


# A = erdos_renyi(100,0.1)
# B = erdos_renyi(100,0.1)

#  edge_list,vertex_attributes = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)

#  SAME.simple_graph_from_edge_list(edge_list)





# edges(er) |> collect


# hypergraph_meta.vertex_attribues



# # #we want to create a list of each clique and each node in each group



# function get_spin_dist(my_vm,my_hypergraph,my_vertex_attribues)
#     bipartite_values = values(my_vertex_attribues)
#     my_N_individuals = count(x->x==1,bipartite_values)
#     my_N_groups = count(x->x==2,bipartite_values)


#     cliques = my_N_individuals+1:my_N_individuals+my_N_groups 

#     nodes_by_clique = map(clique ->neighbors(my_hypergraph,clique),cliques)
#     my_spin_vals = my_vm.spin_network.spin_vals

#     unique_vals =  unique(my_spin_vals)
#     my_count_unique(val,clique_spins) = sum(map(x -> count(x == val),clique_spins))
#     get_values_for_clique(clique_spins_i) = Dict(map(val -> (val,my_count_unique(val,clique_spins_i)),unique_vals))
#     spin_counts = map(clique_i -> get_values_for_clique(my_spin_vals[clique_i]),nodes_by_clique)
#     return spin_counts
# end



begin
    sn = SAME.init_spin_network(graph)
    my_vm = SAME.VoterModel(sn,SAME.voter_model_update)
    vm_ts = SAME.run_model(my_vm,num_steps = 1e4)

    my_hypergraph = hypergraph_meta.hypergraph
    my_vertex_attribues = hypergraph_meta.vertex_attribues

    spin_dist_ts = map( vm_i -> SAME.get_spin_dist(vm_i,my_hypergraph,my_vertex_attribues),vm_ts)
    spin_dist_i = spin_dist_ts[end]
    G_dist = map(spin_dist_i_i -> spin_dist_i_i[1]/sum(values(spin_dist_i_i)),spin_dist_i)
end

histogram(G_dist,bins = 30)
G_dist


# my_sn = vm_ts[1].spin_network
# my_spin_vals = my_sn.spin_vals
# my_graph = my_sn.g
# nodecolordict = Dict(-1=>"lightseagreen", 1 => "orange")
# nodefillc = map(x ->nodecolordict[x],my_spin_vals)
# graph_plot = gplot(my_graph,nodefillc = nodefillc,layout = spring_layout)



# my_binomial(n,u;ϵ = 0.5) = pdf(p_dist,n)*binomial(BigInt(n),BigInt(u))*ϵ^u * (1-ϵ)^(n-u)

# function initialize_u0(N_individuals::Int,p_dist)
#     """
#     initialize_u0 - create an initial distribution for the AME
#     inputs: 
#         N_individuals: the total number of individuals in your populatiojn
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


# t_max = 500
# u₀ = SAME.initialize_u0(N_individuals,p_dist)
# dG = zeros(N_individuals,N_individuals)

# SAME.voter_model_ame!(dG,u₀,0,0)









# p = 0
# t_span = (0.,t_max)
# prob = ODEProblem(SAME.voter_model_ame!,u₀,t_span,p)
# sol = solve(prob,Tsit5(),saveat=1)

# vals = map(x -> sum(x),sol)


# heatmap(sol[6])
# plot!(1:N_individuals,1:N_individuals)
