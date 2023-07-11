using DrWatson
@quickactivate :SAME
using SAME 
# using OrdinaryDiffEq, Plots, RecursiveArrayTools
# using OrdinaryDiffEq:ODESolution
# using GraphPlot
# using Plots
# using PlutoUI


N_individuals = 50
N_p = 0.5
N_groups =  5
M_p = 0.1


begin
    p_dist = Binomial(N_individuals,N_p)                                              
    g_dist = Binomial(N_groups,M_p)
    end


function run_vm_sim(;N_individuals = 50,N_p = 0.5,N_groups = 5,M_p = 0.1,N_steps = 1000,model = "linear"):
    
    #create distributions for p and g
    p_dist = Binomial(N_individuals,N_p)#the hyper-edge size distribution                                              
    g_dist = Binomial(N_groups,M_p)#the hyper-degree distribution
    #initialize the hypergraph
    graph,vertex_attribues = SAME.make_hypergraph(N_groups,N_individuals,p_dist,g_dist) 
	unipartite = SAME.bipartite_projection(graph,vertex_attribues,1)

    sn = SAME.init_spin_network(unipartite)
    my_vm = SAME.VoterModel(sn,SAME.voter_model_update)
    
    N_steps = 10000
    mag_hist = []



end