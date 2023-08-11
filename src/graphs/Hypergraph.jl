
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))#efficent way to unzip a list of tuples



## Data Structures 
#-----------------

abstract type VoterModelGraph end


struct Hypergraph <: VoterModelGraph
    vertices::Set{Int64}
    hyper_edges::Set{Int64}
    edge_list::Dict{Int64, Vector{Int64}}
end

struct BipartiteGraph <: VoterModelGraph
    vertices::Set{Int64}
    edge_list::Vector{Tuple{Int64, Int64}}    
    vertex_attributes::Dict{Int64, Int64}
end


#-----------------
## Construtors
#-----------------

function bipartite_graph(N_groups::Int,N_inds::Int,p_dist,g_dist)
    """
    bipartite_graph(): generates a bipartite_graph in the style of Newman's model in "Community Structure in social and biological networks" 
    inputs:
        N_groups::Int the number of groups or cliques to create
        N_inds::Int the number of individuals to create(may be less than this total)
        p_dist:: The distribution of clique sizes, must be from Distributions.jl
        g_dist:: The distribution of number of cliques belonged to per individual
    output:
        bp_graph::BipartiteGraph- a bipartite graph object
    """
    edge_list,vertex_attributes = generate_hypergraph_bipartite_edge_list(N_groups,N_inds,p_dist,g_dist)
    vertices = keys(vertex_attributes) |> collect
    return BipartiteGraph(Set(vertices),edge_list,vertex_attributes)
end


function hypergraph(N_groups::Int,N_inds::Int,p_dist,g_dist)
    bp_graph=bipartite_graph(N_groups,N_inds,p_dist,g_dist)
    hg = bipartite_to_hypergraph(bp_graph)
    return hg

end


## Simple Graph 
#-----------------
function simple_graph(edge_list::Vector{Tuple{Int64, Int64}})
    """
    simple_graph: make a Graphs.jl SimpleGraph from an edge list
    inputs: 
        edge_list: a list of edges where each edge is a tuple of ints, each int refers to a verticie inde
    returns: 
        graph: a SimpleGraph of the edge_list
    """
    graph = SimpleGraph(length(unique(vcat(unzip(edge_list)...))))

    for edge in edge_list
        add_edge!(graph,edge[1],edge[2])
    end
    return graph
end

function simple_graph(hgraph::Hypergraph)
    """
    simple_graph: make a Graphs.jl SimpleGraph from an edge list
    inputs: 
        hgraph::HyperGraph: a Hypergraph edge list
    returns: 
        graph: a SimpleGraph of the edge_list
    """
    hedge = map(hyperedge -> combinations(hyperedge,2) |> collect,values(hgraph.edge_list))  
    hedge = vcat(hedge...) 
    hedge = map(x->Tuple(x),hedge)
    return simple_graph(hedge)
end

function simple_graph(bp_graph::BipartiteGraph)
    """
    simple_graph: make a Graphs.jl SimpleGraph from an edge list
    inputs: 
        bp_graph::BipartiteGraph: a Bipartite Graph Object 
    returns: 
        graph: a SimpleGraph of the edge_list
    """
    edge_list = bp_graph.edge_list 
    graph = SimpleGraph(length(unique(vcat(unzip(edge_list)...))))

    for edge in edge_list
        add_edge!(graph,edge[1],edge[2])
    end
    return graph
end


function make_simple_graph(N_groups::Int,N_inds::Int,p_dist,g_dist)
    """
    make_hypergraph() - given a set of paramters generate a set a hypergraph consisting of cliques and individuals,
    input: 
        N_groups::Int - the number of groups to generate
        N_inds:: the number of individuals to generate(may be less than this value)
        p_dist::Distributions.jl distribution, the distribution of clique sizes 
        g_dist::Distributions.jl distribution, the distribution of clique sizes
    returns: 
        graph::SimpleGraph, a bipartite graph stored as a unipartite graph with nodes and edges
        my_vertices::Dict{Int64,Int64}, a dictionary with verticie labels and attributes, a vert has an attribute 1 if it is a clique and 2 if it is an individual
    """
    my_edges,vertex_attributes = generate_hypergraph_bipartite_edge_list(N_groups,N_inds,p_dist,g_dist)
    #my_edges = map(x->(Int64(x[1]),Int64(x[2])),my_edges)
    graph = simple_graph(my_edges)
    return graph,vertex_attributes
end


function make_unipartite_community_graph(N_groups::Int,N_inds::Int,p_dist,g_dist)
    # bipartite_graph,vertex_attributes = make_hypergraph(N_groups,N_inds,p_dist,g_dist)#get a bi-partite projection
    edge_list,vertex_attributes = SAME.generate_hypergraph_bipartite_edge_list(N_groups,N_inds,p_dist,g_dist)
    bipartite_graph,vertex_attributes = hypergraph(N_groups,N_inds,p_dist,g_dist) 
    unipartite_graph =  unipartite_projection(bipartite_graph,vertex_attributes,1)

    hypergraph_info = (hypergraph = bipartite_graph,edge_list = edge_list,vertex_attributes = vertex_attributes) 
    return unipartite_graph,hypergraph_info
end


##Generators
##----------

function generate_hypergraph_bipartite_edge_list(N_groups::Int,N_inds::Int,p_dist,g_dist)
    """
    generate_hypergraph_bipartite_edge_list(): generates a hypergraph in the style of Newman's model in "Community Structure in social and biological networks" 
    inputs:
        N_groups::Int the number of groups or cliques to create
        N_inds::Int the number of individuals to create(may be less than this total)
        p_dist:: The distribution of clique sizes, must be from Distributions.jl
        g_dist:: The distribution of number of cliques belonged to per individual
    output:
        edge_list::Vector{Tuple{Int64, Int64}} - the edge list for a bi-partite graph. The first n-indices represent the clique edges and the rest represent individuals
    """

    chairs = []
    butts = []

    #generate chairs
    for i ∈ 1:N_groups
        p_n = rand(p_dist)#select the number of chairs in clique i
        chairs = vcat(chairs,[i for j in 1:p_n])#add p_n chairs belonging to clique i
    end
    


    for i in 1:(N_inds)
            #g_m = rand(g_dist)+1#select the number of butts in clique i
            g_m = rand(g_dist)#select the number of butts in clique i, NOTE may have to fix later if 0s pop up
            g_m = length(butts)+g_m <= length(chairs) ? g_m : length(chairs)-length(butts)#pull a random length or select a length to make the two lists equal if we are bout to go over 
            butts = vcat(butts,[i for j in 1:g_m])#add g_m butts to individuals i
        end
    
    #@infiltrate
    chairs = chairs.+N_inds

    #@infiltrate
    #shuffle the lists    
    chairs_shuff = shuffle(chairs)
    butts_shuff = shuffle(butts)
    #generate edge_list 
    edge_list = zip(chairs_shuff,butts_shuff) |> collect
    edge_list = map(x->(Int64(x[1]),Int64(x[2])),edge_list)
    #create vertex meta_data, if the index is a clique, give it a 0, if the vertex is in individual give it a 1    
    vertex_attributes = Dict( i => i ≤ maximum(butts_shuff) ? 1 : 2 for i ∈ Set(vcat(chairs_shuff,butts_shuff)))
    
    return edge_list,vertex_attributes
    end



##Converters
#----------
function bipartite_to_hypergraph(bp_graph::BipartiteGraph;bipartite_set::Int64 = 1)

    """
    bipartite_to_hypergraph(): converts a bipartite graph into a hypergraph
    inputs:
    bp_graph::BipartiteGraph - a bipartite graph representation of a network you wish to convert to a hypergraph
    bipartite_set::Int - either 0 or 1 the two sets of nodes in the network
    output:
        hgraph::Hypergraph - the converted hypergraph representation of the bipartite graph
    """
    #graph = SAME.simple_graph(bp_graph)
    bipartite_set_vertices = get_bipartite_vertices_type(bp_graph,bipartite_set)
    non_bipartite_vertices = setdiff(Set(bp_graph.vertices),Set(bipartite_set_vertices)) 
    #hyper_edge_list = Dict(hyper_edge => map(x -> x[2], filter(v -> v[1] == hyper_edge, edge_list)) for hyper_edge in bipartite_set_bp_graph.vertices)
    hyper_edge_list = Dict(hyper_edge => map(x -> x[2], filter(v -> v[1] == hyper_edge, bp_graph.edge_list)) for hyper_edge in non_bipartite_vertices)
    #@infiltrate
    hgraph = Hypergraph(bipartite_set_vertices,non_bipartite_vertices,hyper_edge_list)
    return hgraph
end



function unipartite_projection(graph::SimpleGraph, vertex_attributes::Dict{Int64,Int64}, bipartite_set::Int)
    """
    unipartite_projection generate a unipartite projection of a bipartite graph based on one of the bipartite sets.

    Args:
    - graph (SimpleGraph): the input bipartite graph
    - vertex_attributes (Vector{Symbol}): a list of vertex attributes that indicate which of the two bipartite sets a vertex belongs to
    - bipartite_set (Int): the bipartite set to project onto (1 or 2)

    Returns:
    - (SimpleGraph): the unipartite projection of the bipartite graph based on the given bipartite set
    """
    # Find the vertices that belong to the given bipartite set
    my_vertices = collect(vertices(graph))
    bipartite_set_my_vertices = filter(v -> vertex_attributes[v] == bipartite_set, my_vertices)
    non_bipartite_set_my_vertices = setdiff(Set(my_vertices),Set(bipartite_set_my_vertices)) 

    # Generate the unipartite projection of the bipartite graph based on the given bipartite set
    projection = SimpleGraph(length(bipartite_set_my_vertices))


    for vert ∈ non_bipartite_set_my_vertices#for all of the nodes which will be removed in the projection(the cliques)
        bipartite_set_verts = filter(v -> vertex_attributes[v] == bipartite_set, neighbors(graph,vert))#find all the neighbors which are in the unipartite_projection
        new_edges = combinations(bipartite_set_verts,2)#find all the new edges to add 

        map(new_edge ->(add_edge!(projection, new_edge[1],new_edge[2])),new_edges)#add those edges
    end
    return projection
end

function unipartite_projection(hgraph::Hypergraph)
    """
    unipartite_projection() - given a hypergraph construct a unipartite projection: 
        N_groups::Int - the number of groups to generate
        N_inds:: the number of individuals to generate(may be less than this value)
        p_dist::Distributions.jl distribution, the distribution of clique sizes 
        g_dist::Distributions.jl distribution, the distribution of clique sizes
    returns: 
        graph::SimpleGraph, a bipartite graph stored as a unipartite graph with nodes and edges
        my_vertices::Dict{Int64,Int64}, a dictionary with verticie labels and attributes, a vert has an attribute 1 if it is a clique and 2 if it is an individual
    """
    hedge = map(hyperedge -> combinations(hyperedge,2) |> collect,values(hgraph.edge_list))  
    hedge = vcat(hedge...) 
    hedge = map(x->Tuple(x),hedge)
    return SAME.simple_graph(hedge)
end

# function plot_graph(graph;ax =nothing,colors = nothing,layout = nothing)
#     graphplot!(ax,graph,node_color = colors)
# end


#Utilities 
#---------
function get_bipartite_vertices_type(bp_graph::BipartiteGraph,bipartite_set = 1)
    """
    get_bipartite_vertices_type() - returns one of the two sets of verticies in a bipartite graph

    """
    bipartite_set_vertices = Set(filter(v -> bp_graph.vertex_attributes[v] == bipartite_set, bp_graph.vertices))
    return bipartite_set_vertices
end

num_vert(g::Graph)  = nv(g)
num_vert(hg::Hypergraph) = nv(simple_graph(hg))
