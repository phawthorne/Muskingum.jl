using Muskingum
using WatershedSim

baseparams_file = "/Users/hawt0010/Projects/MNRiver/Data/LeSueur/IntegratedModel/data/baseparams.csv"
network_file = "/Users/hawt0010/Projects/MNRiver/Data/LeSueur/IntegratedModel/data/network_table.csv"

wsmodel = StreamModel(baseparams_file, network_file)
X = fill(0.3, wsmodel.nc.n_links)
K = fill(1.0, wsmodel.nc.n_links)
rnmodel = RoutingNetwork(
    wsmodel.nc.outlet_link,
    wsmodel.nc.to_node,
    wsmodel.nc.routing_order,
    X,
    K
)

q_runoff = make_qi(100,100,15)
Δt = 0.01
