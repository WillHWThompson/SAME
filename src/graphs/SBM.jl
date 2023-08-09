
function add_SBM_edges(my_hypergraph::Hypergraph,occupation_prob = 0.1)
    """
    add_sbm_edges(): given a hypergraph add random edges between different cliques, like a stochastic block model with a constant intra-clique occupation rate, these new edges are stored as hyperedges
    inputs:
        my_hypergraph::Hypergraph - a hypergraph you want to add edges to. Note to make it a proper SBM the g_dist must be a dirac delta function with a mean of one
        occupation_prob::Float64 - the occupation probabilty to paramterize the ER graph
    returns 
        ::Hypergraph - a hypergraph with the erdos renyi edges as size two hyper edges
    """
    #create an Erdos-Renyi random graph with the same number of nodes as the hypergraph
    my_unipartite_projection = SAME.simple_graph(my_hypergraph)
    N_inds = length(my_hypergraph.vertices)
    er_graph = erdos_renyi(N_inds,occupation_prob)
    graph_combined = union(my_unipartite_projection,er_graph)
    #extract edge list from ER graph
    er_edge_list = map(x -> (x.src,x.dst),edges(er_graph))
    my_unipartite_projection = SAME.simple_graph(my_hypergraph)
    er_graph = erdos_renyi(N_inds,occupation_prob)

    my_hyperedge_list = my_hypergraph.edge_list
    my_hyperedges = my_hypergraph.hyper_edges
    #add ER edges as new hyperedges and return new ypergraph object 
    hyperedge_counter = maximum(my_hypergraph.hyper_edges)
    @show length(my_hyperedge_list) 
    @show hyperedge_counter
    for edge in edges(er_graph)
        hyperedge_counter+=1
        my_hyperedge_list[hyperedge_counter] = [edge.src,edge.dst]
        push!(my_hyperedges,hyperedge_counter)
    end
    return Hypergraph(my_hypergraph.vertices,my_hyperedges,my_hyperedge_list)
end


function sbm_graph(N_groups,N_inds,p_dist,g_dist,occupation_prob)
    my_hypergraph = hypergraph(N_groups,N_inds,p_dist,g_dist)
    return add_SBM_edges(my_hypergraph,occupation_prob)
end


function print_SBM()
    println("hello world")
end