using Distributed
@everywhere using Distributed
@everywhere println(myid())
@everywhere module Traffic
## imports
using DataStructures
using GeometryBasics
using GraphIO
using Graphs
using MetaGraphsNext
using ProgressMeter
using Statistics
using StatsBase

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
    check_connectedness=false, packets_per_step=0)
    check_connectedness && @show check_connectedness
    #(is_directed(g) ? is_strongly_connected(g) : is_connected(g)) || throw(ArgumentError("Graph must be simply connected"))
    getv = Base.Fix1(getindex, g)
    # to support the order parameter of the continuous case
    active_packets = Int[]
    for t in 1:nsteps
        add_packets!(g, packets_per_step)
        # get queue length for each node, necessary for synchronous update
        pcounts = [length(g[l]) for l in labels(g)]

        push!(active_packets, sum(pcounts))
        # Nodes With Packets
        nwps = findall(!iszero, pcounts)
        # support initial creation mode
        isempty(nwps) && return active_packets, t
        for node in nwps
            # get current packet and remove it from its node
            packet = @lock lck dequeue!(g[node])
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
                pathdist = length(a_star(g, nb, packet))
                queuelen = length(g[nb])
                effdist = distfun(pathdist, queuelen)
                if effdist < mindist
                    mindist = effdist
                    bestnb = nb
                elseif effdist == mindist && rand(Bool)
                    # if there are multiple with the same distance, choose randomly
                    bestnb = nb
                end
            end
            @lock lck enqueue!(g[bestnb], packet)
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
function calculateρ(g, simsteps, avgsteps, τ=1; distfun=(x, y) -> x, packets_per_step=0)
    As = simulate(g, simsteps, distfun; packets_per_step)[1]
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
## simulate for random graphs
#=
g = erdos_renyi(100, 400)
is_connected(g)
mg = initializeSimulationGraph(g, "Erdös-Renyi")
diameter(g) # == 4
add_packets!(mg, 500)
actives, t = simulate(mg)

using GLMakie

ps = collect(0:50:1000)
nmean = 10
Tmaxs = @showprogress map(ps) do p
    samples = map(1:nmean) do _
        add_packets!(mg, p)
        simulate(mg, 1000)
    end
    mean(samples), std(samples)
end

means, stds = first.(Tmaxs), last.(Tmaxs)

plot(ps, means)
errorbars!(ps, means, stds)


dq = simulate(mg, 100000)[2]
npacks = 5
#add_packets!(mg, 5000)
arr = [Point2(0, npacks)]
arrobs = Observable(arr)

plot(arrobs)

nmean = 10
nsteps = 1000
steps = collect(1:nsteps)

collectionarr = zeros(Int, nsteps, nmean)

@showprogress for j = 1:nmean
    # clear left over packets
    simulate(mg, 100000)
    arrobs[] = []
    @showprogress showspeed = true for i = 1:nsteps
        n_active = simulate(mg, 1; packets_per_step=npacks)[1][1]
        push!(arrobs[], (i, n_active))
        notify(arrobs)
        reset_limits!(current_axis())
        yield()
    end
    collectionarr[:, j] = [last(p) for p in arrobs[]]
end

active_mean = mean(collectionarr, dims=2)[:, 1]
active_std = std(collectionarr, dims=2)[:, 1]
As = map(t -> -(t...), zip(active_mean[end-100:end], active_mean[end-200:end-100]))
mean(As)

fig, ax, plt = errorbars(steps, active_mean, active_std, color=:orange)
plot!(steps, active_mean)
ax.xlabel = "Time t"
ax.ylabel = "Number of active packets A"
ax.title = "Traffic model on a single instance of an ER graph with 100 nodes and 400 edges"
save("ER_100_400_3pps.png", fig, px_per_unit=3)
=#
## real data AS map
@everywhere begin
    asg = Graph(loadgraph("as19981229.txt", GraphIO.EdgeList.EdgeListFormat()))
    simg = initializeSimulationGraph(asg, "AS level internat map, 29/12/1998")
    is_connected(simg)

    ppss = collect(1:1:15)
    hs = 0.7:0.1:1
    npasses = 10
    simsteps = 1000
    avgsteps = 50
    τ = 1
end
fig = Figure()
ax = Axis(fig[1, 1])

@everywhere function updatefig(ρs, Δs, h)
    @info "calling updatefig on $(myid())"
    plot!(ppss, ρs, label="h=$h")
    errorbars!(ppss, ρs, Δs, label="h=$h")
    autolimits!(ax)
end

display(fig)

println("starting loop")
@sync @distributed for h in hs
    distf = getEffDist(h)
    ρs = Float64[]
    Δs = Float64[]
    locg = initializeSimulationGraph(asg, "AS level internat map, 29/08/1999")
    @info "starting for h=$h"
    @showprogress "PPS loop" for pps in ppss
        @info "currently at $pps for h=$h"
        Σρ = 0.0±0.0
        lck = ReentrantLock()
        Threads.@threads for pass in 1:npasses
            @info "$pass/$npasses for $pps and $h on process $(myid())T$(Threads.threadid())"
            @time ρ, Δρ = calculateρ(locg, simsteps, avgsteps, τ; distfun=distf, packets_per_step=pps)
            @lock lck Σρ += ρ±Δρ
        end
        Σρ /= npasses
        ρ = Σρ.val
        Δρ = Σρ.err
        push!(ρs, ρ)
        push!(Δs, Δρ)
    end
    @spawnat 1 updatefig(ρs, Δs, h)
end

ax.title = "Order Parameter ρ for multiple routing strategies"
ax.xlabel = "Packets added per step"
ax.ylabel = "Order parameter ρ"
fig[1, 2] = Legend(fig, ax, merge=true)

save("reduced_plot_h_pps.png", fig)


## use real data, padua street network
# see `traffic_padova.jl`