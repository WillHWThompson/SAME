function make_SBM(N_groups,N_individuals,p_dist,g_dist,occupation_prob)
    edge_list,vertex_attribues = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_individuals,p_dist,g_dist)
    hypergraph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist) 

	unipartite = SAME.bipartite_projection(hypergraph,vertex_attribues,2)

    ER_graph = erdos_renyi(N_individuals,occupation_prob)
    graph_combined = union(unipartite,ER_graph)
    hypergraph_info = (hypergraph = hypergraph,edge_list = edge_list,vertex_attribues = vertex_attribues) 
    return graph_combined,hypergraph_info
end

function print_SBM()
    println("hello world")
end