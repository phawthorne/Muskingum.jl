using Printf
using Plots
using Distributions

function simple_example()
    # make input timeseries
    d = Gamma(4, 1)
    T = 30
    L = 5
    N = 2000
    M = 100
    t = range(0, T, length=N)

    Δt = step(t)
    X = 0.3
    K = 1.0

    qi = fill(0.0, (L, N))
    qe = fill(0.0, (L, N))
    qi[1,:] = 100 .* pdf(d, t)

    c0 = (-2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    c1 = (2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    c2 = (2*(1-X) - Δt/K) / (2*(1-X) + Δt/K)

    for i in 1:(N-1)
        qe[1, i+1] = c0 * qi[1, i+1] + c1 * qi[1, i] + c2 * qe[1, i]
        for l in 2:L
            qe[l, i+1] = c0 * qe[l-1, i+1] + c1 * qe[l-1, i] + c2 * qe[l, i]
        end
    end

    plot(t, qi[1,:], label="input")
    for l in 1:L
        plot!(t, qe[l, :], label=@sprintf "output %d" l)
    end
    plot!()
end


function reach_params_example()
    d = Gamma(4, 1) # input distribution
    T = 30          # duration
    L = 5           # number of links
    N = 2000        # number of timesteps
    M = 100         # input scale param

    t = range(0, T, length=N)
    Δt = step(t)

    X = [0.25, 0.3, 0.1, 0.4, 0.3]
    K = [2.0, 1.8, 1.5, 1.0, 1.5]

    qi = fill(0.0, (L, N))
    qe = fill(0.0, (L, N))
    qi[1,:] = 100 .* pdf(d, t) + 80 .* pdf(d, t .- 12)

    c0 = fill(0.0, L)
    c1 = fill(0.0, L)
    c2 = fill(0.0, L)
    @. c0 = (-2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    @. c1 = (2*X + Δt/K) / ( 2*(1-X) + Δt/K)
    @. c2 = (2*(1-X) - Δt/K) / (2*(1-X) + Δt/K)

    for i in 1:(N-1)
        qe[1, i+1] = c0[1] * qi[1, i+1] + c1[1] * qi[1, i] + c2[1] * qe[1, i]
        for l in 2:L
            qe[l, i+1] = c0[l] * qe[l-1, i+1] + c1[l] * qe[l-1, i] + c2[l] * qe[l, i]
        end
    end

    plot(t, qi[1,:], label="Qi1")
    for l in 1:L
        plot!(t, qe[l, :], label=@sprintf "Qe%d" l)
    end
    title!("Muskingum")
    xlabel!("Time")
    ylabel!("Flow")
    plot!()
end


function single_path_sim(T, nsteps, qi1::Vector{Float64},
                         X::Vector{Float64}, K::Vector{Float64})
     nlinks = length(X)
     t = range(0, T, length=nsteps)
     Δt = step(t)

     qi = fill(0.0, (nlinks, nsteps))
     qe = fill(0.0, (nlinks, nsteps))
     qi[1, :] = qi1

     c0 = fill(0.0, nlinks)
     c1 = fill(0.0, nlinks)
     c2 = fill(0.0, nlinks)
     @. c0 = (-2*X + Δt/K) / ( 2*(1-X) + Δt/K)
     @. c1 = (2*X + Δt/K) / ( 2*(1-X) + Δt/K)
     @. c2 = (2*(1-X) - Δt/K) / (2*(1-X) + Δt/K)

     for i in 1:(nsteps-1)
         qe[1, i+1] = c0[1] * qi[1, i+1] + c1[1] * qi[1, i] + c2[1] * qe[1, i]
         for l in 2:nlinks
             qe[l, i+1] = c0[l] * qe[l-1, i+1] + c1[l] * qe[l-1, i] + c2[l] * qe[l, i]
         end
     end

     return qi, qe
 end


function make_qi(days, steps_per_day, event_spacing)

    M = Uniform(50, 200)
    d = Gamma(3, 1) # input distribution
    S = Poisson(event_spacing)

    qi = fill(0.0, days*steps_per_day)
    T = 1
    te = range(0, 10, length=10*steps_per_day)
    while T + 10*steps_per_day < days*steps_per_day
        event_start = T
        event_end = T + 10*steps_per_day-1
        event_q = rand(M) .* pdf(d, te)
        qi[event_start:event_end] .+= event_q
        T +=  rand(S) * steps_per_day
    end

    return qi
end
