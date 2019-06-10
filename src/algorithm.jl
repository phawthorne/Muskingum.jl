using Distributions
using Parameters
using Plots
using Printf

@with_kw struct RoutingNetwork
    outlet_link::Int64
    to_node::Array{Int64,1}
    routing_order::Array{Int64,1}
    X::Array{Float64,1}
    K::Array{Float64,1}
end

@with_kw struct RoutingNetwork2
    outlet_link::Int64
    to_node::Array{Int64,1}
    routing_order::Array{Int64,1}
    X::Array{Float64,1}
    K::Array{Float64,1}
    qi::Array{Float64,2}
    qe::Array{Float64,2}
end

function muskingum_routing(network::RoutingNetwork, q_runoff::Array{Float64, 2},
                           Δt::Float64)
    @unpack to_node, routing_order, X, K = network

    n_links = length(to_node)
    T = size(q_runoff, 2)

    qi = deepcopy(q_runoff)
    qe = fill(0.0, size(q_runoff))
    qe[:,1] = qi[:,1]

    c0 = fill(0.0, n_links)
    c1 = fill(0.0, n_links)
    c2 = fill(0.0, n_links)
    @. c0 = (-2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    @. c1 = (2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    @. c2 = (2*(1-X) - Δt/K) / (2*(1-X) + Δt/K)

    for t in 1:(T-1)
        for l in routing_order
            qe[l, t+1] = c0[l]*qi[l,t+1] + c1[l]*qi[l,t] + c2[l]*qe[l,t]
            if to_node[l] >= 0
                qi[to_node[l], t+1] += qe[l, t+1]
            end
        end
    end

    return qi, qe
end


function muskingum_routing!(network::RoutingNetwork2,
                            q_runoff::Array{Float64,2},
                            Δt::Float64)
    @unpack to_node, routing_order, X, K, qi, qe = network

    n_links = length(to_node)
    T = size(q_runoff, 2)

    # reset qi, qe
    @. qi = q_runoff
    @. qe = 0.0

    c0 = fill(0.0, n_links)
    c1 = fill(0.0, n_links)
    c2 = fill(0.0, n_links)
    @. c0 = (-2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    @. c1 = (2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    @. c2 = (2*(1-X) - Δt/K) / (2*(1-X) + Δt/K)

    for t in 1:(T-1)
        for l in routing_order
            qe[l, t+1] = c0[l]*qi[l,t+1] + c1[l]*qi[l,t] + c2[l]*qe[l,t]
            if to_node[l] >= 0
                qi[to_node[l], t+1] += qe[l, t+1]
            end
        end
    end

end
