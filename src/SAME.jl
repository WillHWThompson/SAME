module SAME
using Graphs
using Reexport
@reexport using Random,Distributions,Graphs,Combinatorics
@reexport using Revise

export models,utils,graphs

 include("models/BaseModels.jl")
 include("models/IsingModel.jl")
 include("models/VoterModel.jl")

 println("Including Graphs")
 include("graphs/Hypergraph.jl")
 include("graphs/SBM.jl")

 include("utils/ModelUtils.jl")
 println("All files included")
 include("ame/AME.jl")

 greet() = println("hello world2")


end # module SAME
