using DrWatson
@quickactivate :SAME
using SAME 
using GraphPlot
N_individuals = 500
N_groups = 10  


#p_dist = Binomial(N_individuals,0.1)                                              
#g_dist = Binomial(N_groups,0.5)

p_dist = Dirac(N_individuals)                                              
g_dist = Dirac(N_groups)


edge_list,vertex_attribues = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)
graph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist) 

sn = SAME.init_spin_network(graph)
my_vm = SAME.VoterModel(sn,SAME.voter_model_update)
vm_ts = SAME.run_model(my_vm)


my_sn = vm_ts[1].spin_network
my_spin_vals = my_sn.spin_vals
my_graph = my_sn.g
nodecolordict = Dict(-1=>"lightseagreen", 1 => "orange")
nodefillc = map(x ->nodecolordict[x],my_spin_vals)
graph_plot = gplot(my_graph,nodefillc = nodefillc,layout = circular_layout)

