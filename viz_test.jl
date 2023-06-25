using DrWatson, Revise
@quickactivate :SAME
using .SAME

N_individuals = 500
N_groups = 10

p_dist = Binomial(N_groups,0.01)
g_dist = Binomial(N_individuals,0.05)


rand(p_dist)
edge_list,vertex_attribues = generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)
graph,vertex_attribues = make_hypergraph(N_groups,N_individuals,p_dist,g_dist)

