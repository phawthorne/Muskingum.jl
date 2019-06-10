using Muskingum

qi_onedim = make_qi(11, 20, 15)

linear_network = RoutingNetwork(
    5,
    [2,3,4,5,-1],
    [1,2,3,4,5],
    [0.2, 0.2, 0.2, 0.2, 0.2],
    [1.0, 1.0, 1.0, 1.0, 1.0],
)

qi_sp, qe_sp = single_path_sim(11, length(qi_onedim), qi_onedim,
                               linear_network.X,
                               linear_network.K)

qi_mat = fill(0.0, (5, length(qi_onedim)))
qi_mat[1,:] = qi_onedim
qi_net, qe_net = network_muskingum(linear_network, qi_mat, 0.050228)


branched_network = RoutingNetwork(
    7,
    [5,5,6,6,7,7,-1],
    [1,2,3,4,5,6,7],
    fill(0.2, 7),
    fill(1.0, 7),
)

qi_branched = fill(0.0, (7, length(qi_onedim)))
for i in 1:4
    qi_branched[i, :] .= qi_onedim
end
qi_b, qe_b = muskingum_routing(branched_network, qi_branched, 0.050228)

bn2 = RoutingNetwork2(
    7,
    [5,5,6,6,7,7,-1],
    [1,2,3,4,5,6,7],
    fill(0.2, 7),
    fill(1.0, 7),
    fill(0.0, (7, length(qi_onedim))),
    fill(0.0, (7, length(qi_onedim)))
)

muskingum_routing!(bn2, qi_branched, 0.050228)
