using Distributed
@everywhere using Distributed
@everywhere println(myid())
@everywhere module Traffic
## imports
using DataStructures
using GeometryBasics
using GraphIO
using Graphs
using Graphs: weights
using MetaGraphsNext
using ProgressMeter
using Statistics
using StatsBase
using SparseArrays

export initializeSimulationGraph, simulate, calculateρ, getEffDist

## change directory
isdir("code") && cd("code/task_16/")

## simulation functions

"""
    initializeSimulationGraph(g, title="")

initializes a `MetaGraph` as expected by the simulation from a simple `Graph` g
with `graph_data` set to `title`
"""
function initializeSimulationGraph(g, title="")
    vs = vertices(g)
    vertices_description = [v => Queue{Int}() for v in vertices(g)]
    edges_description = [(src(e), dst(e)) => nothing for e in edges(g)]
    MetaGraph(g, vertices_description,
        edges_description, title, _ -> 1.0, 0.0
    )
end  # function initializeSimulationGraph

"""
    getEffDist(h=1.0)
    Return the effective distance function for a given `h`.
"""
function getEffDist(h=1.0)
    h > 1 && error("h needs to be smaller or equal to 1 but is $h")
    if isone(h)
        (x, _) -> x
    else
        (x, y) -> h * x + (1 - h) * y
    end
end  # function getEffDist

"""
    add_packets!(g, n=100)
    adds `n` packets to the graph `g`
"""
function add_packets!(g, n=100)
    for _ in 1:n
        # choose random source
        s = d = rand(vertices(g))
        # ensure destination != source
        while d == s
            d = rand(vertices(g))
        end
        enqueue!(g[s], d)
    end
    g
end  # function add_packets!

"""
    simulate(g, nsteps=100, distfun=(x, y)->x;
        check_connectedness=false, packets_per_step=0)

Simulate a simple traffic simulation on graph `g` for `nsteps` with a changeable function `distfun` used for
    optimal route calculation.
    
# Arguments:
- `g`: Graph, needs to be (strongly) connected, each vertex holds a queue of packets
- `nsteps`: timesteps to simulate the packet progression for
- `distfun`: effective distance function, takes shortest path length and queue length
    default: just shortest path
- `check_connectedness`: if function should check the (strong) connectedness of the graph
- `packets_per_step`: amount of packets to add to the graph at each timesteps

# Notes
can be used both in first fill, then simulate mode and in continuous instream mode
"""
function simulate(g, nsteps=100, distfun=getEffDist(1);
    check_connectedness=false, packets_per_step=0, distfun_elemwise=true)
    check_connectedness && @show check_connectedness
    #(is_directed(g) ? is_strongly_connected(g) : is_connected(g)) || throw(ArgumentError("Graph must be simply connected"))
    getv = Base.Fix1(getindex, g)
    # to support the order parameter of the continuous case
    active_packets = Int[]

    orig_distmat = sparse(weights(g))

    for t in 1:nsteps
        add_packets!(g, packets_per_step)
        # get queue length for each node, necessary for synchronous update
        pcounts = [length(g[l]) for l in labels(g)]

        push!(active_packets, sum(pcounts))
        # Nodes With Packets
        nwps = findall(!iszero, pcounts)
        # support initial creation mode
        isempty(nwps) && return active_packets, t

        distmat = copy(orig_distmat)
        if distfun_elemwise
            for i in axes(distmat, 1)
                for j in axes(distmat, 2)
                    dij = distmat[i, j]
                    if !iszero(dij)
                        distmat[i, j] = distfun(dij, pcounts[j])
                    end
                end
            end
        end

        for node in nwps
            # get current packet and remove it from its node
            packet = dequeue!(g[node])
            # neighbours
            nbs = outneighbors(g, node)
            # the destination is in the neighbourhood
            if any(==(packet), nbs)
                # -> remove the packet from circulation by not inserting it anywhere
                continue
            end
            mindist = Inf
            bestnb = rand(nbs)
            for nb in nbs
                pathdist = length(a_star(g, nb, packet, distmat))
                queuelen = length(g[nb])
                # if all the individual weights are adjusted,
                # take this effective distance unchanged
                # conditional only used for figure in the main part, not the appendix
                effdist = distfun_elemwise ? pathdist : distfun(pathdist, queuelen)
                if effdist < mindist
                    mindist = effdist
                    bestnb = nb
                elseif effdist == mindist && rand(Bool)
                    # if there are multiple with the same distance, choose randomly
                    bestnb = nb
                end
            end
            enqueue!(g[bestnb], packet)
        end
    end
    return active_packets, nsteps
end  # function simulate

"""
    calculateρ(g, simsteps, avgsteps, τ=1; distfun=(x, y)->x, packets_per_step=0)
    Calculate the order parameter ρ for a given graph `g` and simulation parameters.

    # Arguments:
    - `simsteps`: # of steps to run the simulation for
    - `avgsteps`: last `avgsteps` of the simulation get used to calculate `ρ`
    - `τ`: time to wait to calculate a single `ρ` samples

    See also for the rest [`simulate`](@ref)
"""
function calculateρ(g, simsteps, avgsteps, τ=1, distfun_elemwise=true; distfun=(x, y) -> x, packets_per_step=0)
    As = simulate(g, simsteps, distfun; packets_per_step, distfun_elemwise)[1]
    r1 = lastindex(As)-avgsteps:lastindex(As)
    r2 = r1 .- τ
    ρs = float(map((x, y) -> x - y, As[r1], As[r2]))
    ρs ./= τ * packets_per_step
    avg = mean(ρs)
    stdev = std(ρs; mean=avg)
    avg, stdev
end  # function calculateρ
end #Module
@everywhere using .Traffic, GraphIO, Graphs, GLMakie, ProgressMeter, Measurements

## real data AS map
@everywhere begin
    asg = DiGraph(Graph(loadgraph("as19990829.txt", GraphIO.EdgeList.EdgeListFormat())))
    simg = initializeSimulationGraph(asg, "AS level internat map, 29/08/1999")
    is_connected(simg)

    ppss = collect(1:1:15)
    hs = [0.5, 0.5, 0.7, 0.7, 0.9, 0.9, 1]
    elemwise = [true, false, true, false, true, false, false]
    npasses = 10
    simsteps = 500
    avgsteps = 50
    τ = 1
end
fig = Figure(size=(1000, 500))
ax = Axis(fig[1, 1])

@everywhere function updatefig(ρs, Δs, h, elem)
    @info "calling updatefig on $(myid())"
    plot!(ppss, ρs, label="h=$h, elemwise=$elem")
    errorbars!(ppss, ρs, Δs, label="h=$h, elemwise=$elem")
    autolimits!(ax)
end

display(fig)

println("starting loop")
@sync @distributed for (h, elem) in collect(zip(hs, elemwise))
    distf = getEffDist(h)
    ρs = Float64[]
    Δs = Float64[]
    @info "starting for h=$h"
    @showprogress "PPS loop" for pps in ppss
        @info "currently at $pps for h=$h"
        Σρ = 0.0±0.0
        lck = ReentrantLock()
        Threads.@threads for pass in 1:npasses
            locg = initializeSimulationGraph(asg, "AS level internat map, 29/12/1998")
            @info "$pass/$npasses for $pps and $h on process $(myid())T$(Threads.threadid())"
            @time ρ, Δρ = calculateρ(locg, simsteps, avgsteps, τ, elem; distfun=distf, packets_per_step=pps)
            @lock lck Σρ += ρ±Δρ
        end
        Σρ /= npasses
        ρ = Σρ.val
        Δρ = Σρ.err
        push!(ρs, ρ)
        push!(Δs, Δρ)
    end
    @spawnat 1 updatefig(ρs, Δs, h, elem)
end

ax.title = "Order Parameter ρ for multiple routing strategies"
ax.xlabel = "Packets added per step"
ax.ylabel = "Order parameter ρ"
fig[1, 2] = Legend(fig, ax, merge=true)

save("reduced_plot_h_elemwise_pps_notfirst.png", fig)


## use real data, padua street network
# see `traffic_padova.jl`