module Muskingum

include("algorithm.jl")
export RoutingNetwork, RoutingNetwork2, muskingum_routing, muskingum_routing!

include("Tmp.jl")
export simple_example, reach_params_example, single_path_sim,
       make_qi

end # module
