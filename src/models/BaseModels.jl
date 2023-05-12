
using DrWatson
import Base.@kwdef
@quickactivate "SpinModelsAME"


abstract type  SpinModel end

mutable struct SpinNetwork
    """
    SpinNetwork a datastrcture that contains a network of spins
    inputs:
        spin_vals: a vector with length equal to the number of nodes in the network assigning a spin to each node
        g: a graph structure that is used in the network
        spin_meta_data: a vector of dictionaries containing any spin specific meta-data
    """
    spin_vals::Vector{Int64}
    g::Graph
    spin_meta_data::Vector{Any}
end


function init_spin_network(network;spin_pmf_dict= Dict(-1 =>1/2, 1 => 1/2),spin_meta_data = Vector())
    """
    init_spin_network: a constructor the creates a SpinNetwork struct

    inputs: 
        network: a Graphs.jl network the represetns the topology of the spin network 
        spin_pmf_dict: a dictioanry with keys as the support for a pmf and values and the probabilites of that support, the values represent the posible spin states
        spin_meta_data: a vector of dictionaries each containing meta_data about the spins of each node
    """
    xs,ps = collect(keys(spin_pmf_dict)),collect(values(spin_pmf_dict))#extract support and probs from dictionary
    pmf = DiscreteNonParametric(xs,ps)#init nonparametic pmf to draw from

    spin_vals = rand(pmf,nv(network))#assign each node in the network a spin, with a prob drawn from the pmf you defined

    sn = SpinNetwork(spin_vals,network,spin_meta_data)
    return sn
end





function get_spin_dist(my_vm,my_hypergraph,my_vertex_attribues)
    bipartite_values = values(my_vertex_attribues)
    my_N_individuals = count(x->x==1,bipartite_values)
    my_N_groups = count(x->x==2,bipartite_values)


    cliques = my_N_individuals+1:my_N_individuals+my_N_groups-3

    nodes_by_clique = map(clique ->neighbors(my_hypergraph,clique),cliques)
    my_spin_vals = my_vm.spin_network.spin_vals

    unique_vals =  unique(my_spin_vals)
    my_count_unique(val,clique_spins) = sum(map(x -> count(x == val),clique_spins))
    get_values_for_clique(clique_spins_i) = Dict(map(val -> (val,my_count_unique(val,clique_spins_i)),unique_vals))
    spin_counts = map(clique_i -> get_values_for_clique(my_spin_vals[clique_i]),nodes_by_clique)
    return spin_counts
end


function run_model(my_model::SpinModel;num_steps = 1000)
    println(num_steps)
    history = map(x ->my_model.update_rule(my_model),1:num_steps)
    return history
end