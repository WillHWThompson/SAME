
using DrWatson
""""
----------------------
ISING MODEL FUNCTIONS
---------------------
"""

@kwdef mutable struct IsingModel  <: SpinModel
    spin_network::SpinNetwork
    update_rule::Function
    β::Float64#inverse temperature 
    h::Float64#external magnetic field
    J::Float64#the strength of the coupling

end

function metropolis_hastings_update(my_ising_model::IsingModel)
    """
    An implementation of the Metropolis-Hastings algorithm for simulating the single spin flip dynamics of the ising model for a single spin flip
    input: 
        my_ising_model::IsingModel
    output:
        my_ising_model::Ising Model
    refs: 
        https://en.wikipedia.org/wiki/Ising_model#Monte_Carlo_methods_for_numerical_simulation
    """
    my_ising_model = deepcopy(my_ising_model)
    spin_idx = rand(1:nv(my_ising_model.spin_network.g))
    spin = my_ising_model.spin_network.spin_vals[spin_idx]

    # calculate the energy of the spin flip
    energy_diff = 2*my_ising_model.J*spin*sum(my_ising_model.spin_network.spin_vals[neighbors(my_ising_model.spin_network.g, spin_idx)])
    energy_diff += -2 * my_ising_model.h * spin

    #decide to accept or reject the change
    if energy_diff <= 0 || rand() < exp(-energy_diff / my_ising_model.β)
        # accept the spin flip
        my_ising_model.spin_network.spin_vals[spin_idx] = -spin
    end
    return my_ising_model
end


function get_Hamiltonian(my_ising_model::IsingModel)
    """
    get_Hamiltonian: calculate the Hamiltonian for the current configuration of the Ising Model 
    input: 
        my_ising_model::IsingModel
    output: 
        H::Float64 - the value of the Hamiltonian. 
    """
    H = 0.0
    for i in eachindex(my_ising_model.spin_network.spin_vals)
        h = my_ising_model.h
        J = my_ising_model.J
        
        spin = my_ising_model.spin_network.spin_vals[i]
        H-= h*spin
        for neighbor_idx in neighbors(my_ising_model.spin_network.g,i)
            H -= J*spin*my_ising_model.spin_network.spin_vals[neighbor_idx] 
        end
    end
    return H 
end

function get_magnetization(my_ising_model::IsingModel)
    """
    get_Magnetization: returns the magnetization of the current configuration of the ising model
    input: 
        my_ising_model::IsingModel
    output: 
        M::Float64 - the magnetization
    """
    M = sum(my_ising_model.spin_network.spin_vals)/length(my_ising_model.spin_network.spin_vals)
    return M
end

function IsingModel_1_step(my_ising_model::IsingModel)
    new_ising_model = deepcopy(metropolis_hastings_update(my_ising_model))
    # new_h = get_Hamiltonian(new_ising_model)
    # new_m = get_magnetization(new_ising_model)
   #vals = Dict(:ising_model => new_ising_model, :h => new_h, :m => new_m)
    return new_ising_model 
end


function plot_spin_network(spin_network::SpinNetwork)
    x = [i[1] for i in spin_network.sites]
    y = [i[2] for i in spin_network.sites]
    z = [spin_network[i] for i in 1:length(spin_network)]

    scatter(x, y, zcolor=z, marker=(:square, 10), legend=false, colorbar=false)
    plot!(xlims=(minimum(x)-0.5, maximum(x)+0.5), ylims=(minimum(y)-0.5, maximum(y)+0.5), aspect_ratio=1)
end