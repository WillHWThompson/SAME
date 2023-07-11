
using DrWatson
""""
----------------------
ISING MODEL FUNCTIONS
---------------------
"""

function test2()
    println("test")
end


abstract type VoterModel <: SpinModel end

#mutable struct LinearVoterModel{S <: SpinNetwork}  <: SpinModel
mutable struct LinearVoterModel{S <: SpinNetwork} 
"""
VoterModel: A strcut to store a voter model 
inputs: 
    spin_network::SpinNetwork - a spin network topology with inital values and a graph
    update_rule::Function - the rule used to update the graph
"""
    spin_network::S
    update_rule::Function
end

""""
----------------------
Constructors
---------------------
"""

function linear_voter_model(spin_model::SpinNetwork)
    return LinearVoterModel(spin_model,linear_voter_model)
end



function linear_voter_model_update(my_voter_model::VoterModel)
    """
    An implementation of the voter model update rule, each generation a single new spin is chosen to update
    input: 
        my_voter_model::VoterModel
    output:
        my_voter_model::VoterModel
    refs: 
        http://physics.bu.edu/~redner/542/book/spin.pdf
    """
    my_voter_model = deepcopy(my_voter_model)
    spin_idx = rand(1:nv(my_voter_model.spin_network.g))

    neighbor_spins = my_voter_model.spin_network.spin_vals[neighbors(my_voter_model.spin_network.g, spin_idx)]#get the neighboring spins
   
    n_neighbors = length(neighbor_spins)
    if n_neighbors > 0
        new_spin = neighbor_spins[rand(1:n_neighbors)]#select a random spin value from neighbors
        my_voter_model.spin_network.spin_vals[spin_idx] = new_spin#select spin
    end
    return my_voter_model
end


function get_magnetization(my_voter_model::VoterModel)
    """
    get_Magnetization: returns the magnetization of the current configuration of the voter model
    input: 
        my_voter_model::voterModel
    output: 
        M::Float64 - the magnetization
    """
    M = sum(my_voter_model.spin_network.spin_vals)/length(my_voter_model.spin_network.spin_vals)
    return M
end


function run_model(d::Dict)
    @unpack N_inds, N_p, N_groups, M_p, N_steps = d
    p_dist = Binomial(N_inds,N_p)                                              
    g_dist = Binomial(N_groups,M_p)
    hg = hypergraph(N_groups,N_inds,p_dist,g_dist)
    sn = init_spin_network(hg)
    vm = linear_voter_model(sn)
    for i in 1:N_steps
        vm.update_rule(vm.spin_network)
    end
    calculate_hyperedge_spin_dist(sn)
end