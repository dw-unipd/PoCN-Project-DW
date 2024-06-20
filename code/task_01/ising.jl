## imports
using Graphs
using Statistics
using ProgressMeter

## function definitions

include("../common_utils.jl")

"""
    Hamilton(g, statevec; J=1.0, H=0.0)
    Find the energy of a state of a graph.
    g: a graph
    statevec: vector of spin states ±1, indexed by node id
"""
function Hamilton(g, statevec; J=1.0, H=0.0)
   nv(g) == length(statevec) || error("node count of graph and length of state vector don't match. $(nv(g)) != $(length(statevec))")
   E = 0.0
   @inbounds for e in edges(g)
      s, d = src(e), dst(e)
      E += -J * statevec[s] * statevec[d]
   end
   iszero(H) && return E
   @inbounds for n in vertices(g)
      E += -H * statevec[n]
   end
   E
end  # function Hamilton


"""
   magnetization(statevec)
   calculates average absolute magnetization
"""
magnetization(statevec) = sum(statevec)

"""
   sweep!(g, statevec, T)
   on average tests for each spin, runs metropolis sample
   g: Graph
   statevec: vector of spins
   T: current temperature
"""
function sweep!(g, statevec, T; probf=expProb)
   nrange = 1:nv(g)
   for i in nrange
      n = rand(nrange)
      s = statevec[n]
      dE = 0
      # ΔE = -2 * Σ_j s_i*s_j
      @inbounds for neigh in neighbors(g, n)
         dE += statevec[neigh]
      end
      dE *= 2 * s
      if rand() < probf(dE, T)
         statevec[n] = -s
      end
   end
end  # function sweep  

"""
   simulated_annealing(g, statevec;
   initial_temp=1.0e1, max_steps=10 * nv(g), probf=expProb,
   adjust_temp=true)
   Find a Hamiltonian state using simulated annealing.
   g: a graph
   statevec: vector of spin states ±1, indexed by node id
   initial_temp: initial temperature
   max_steps: maximum number of steps
   probf: function to compute probability of accepting a new state
"""
function simulated_annealing(g, statevec;
   initial_temp=1.0e1, end_temp=0.0, max_steps=100, probf=expProb,
   adjust_temp=true)
   if initial_temp < end_temp
      @warn "setting initial_temp to end_temp adn disabling temperature scaling 
      because initial_temp $initial_temp < $end_temp end_temp"
      initial_temp = end_temp
      adjust_temp = false
   end
   T = initial_temp
   end_temp = iszero(end_temp) ? initial_temp * 1e-10 : end_temp
   Tdamp = exp(log(end_temp / initial_temp) / max_steps)
   E = Hamilton(g, statevec)
   # accumulator for E, E^2, M, M^2
   # helps with auxiliary 
   accum = zeros(4)
   @inbounds for i in 1:max_steps
      # linear temperature descent
      adjust_temp && (T *= Tdamp)
      sweep!(g, statevec, T)
      E = Hamilton(g, statevec)
      M = magnetization(statevec)
      accum[1] += E
      accum[2] += E^2
      accum[3] += (M)
      accum[4] += M^2
   end
   #println(accum)
   # normalization by amount of steps
   accum ./= max_steps
   # also needs to be normalized by number of contributers, square values doubly so
   accum[1:2] ./= ne(g)
   accum[2] /= ne(g)
   accum[3:4] ./= nv(g)
   accum[4] /= nv(g)
   return statevec, accum
end  # function simulated_annealing

function equilibrate(g, statevec, T=1.0, steps=100)
   simulated_annealing(g, statevec;
      initial_temp=T, max_steps=steps, adjust_temp=false)
end  # function equilibrate

function expProb(dE, T)
   return exp(-dE / T)
end  # function expProb

function neighbour!(x, x_proposed)
   copy!(x_proposed, x)
   n = n = rand(1:length(x))
   x_proposed[n] = -x[n]
   nothing
end  # function neighbour!

## plotting basic behaviour

using GraphMakie, CairoMakie
#=
sz = 500
sbm = scale_free_with_cutoff(sz, 3.2, 2) 
sv = rand(-1.0:2:1, nv(sbm))
preE = Hamilton(sbm, sv)
graphplot(sbm, node_color=sv, edge_width=0.1)
show(current_figure)
simulated_annealing(sbm, sv; max_steps=1e3, initial_temp=1.0, adjust_temp=false)
postE = Hamilton(sbm, sv)
graphplot(sbm, node_color=sv, edge_width=0.1)
=#

sz = 1000
nsamples = 100
Ts = collect(0:0.1:5)

Es = [zeros(nsamples) for _ = Ts]
mags = [zeros(nsamples) for _ = Ts]
χs = [zeros(nsamples) for _ = Ts]
Cvs = [zeros(nsamples) for _ = Ts]


@showprogress for (id, T) in enumerate(Ts)
   Threads.@threads for s in 1:nsamples
      g = static_scale_free(sz, 5 * sz, 2.0)
      sv = ones(nv(g))
      # equilibration for 50 iterations
      simulated_annealing(g, sv; initial_temp=T + 1, end_temp=T, max_steps=100)
      # now with extraction of thermodynamic quantities
      sv, accum = equilibrate(g, sv, T, 100)
      Es[id][s] = accum[1]
      mags[id][s] = accum[3]
      # χ = Var(M)/T
      χs[id][s] = (accum[4] - accum[3]^2) / T
      # C_v = Var(E)/T^2
      Cvs[id][s] = (accum[2] - accum[1]^2) / T^2
   end
   GC.gc()
end
Ē = abs.(mean.(Es))
varE = std.(Es)
plot(Ts, Ē, markersize=10, label="E")
errorbars!(Ts, Ē, varE, whiskerwidth=5)

maḡ = abs.(mean.(mags))
varmag = std.(mags)
plot!(Ts, maḡ, markersize=10, label="M")
errorbars!(Ts, maḡ, varmag, whiskerwidth=5)

current_axis().title = "Ising Model behaviour on scale free network with PL exponent 2.0"
current_axis().xlabel = "Temperature T"
axislegend()
save("Ising_EM_SF_20.pdf", current_figure())

χ̄ = mean.(χs)
χ̄ ./= maximum(χ̄)
varχ = std.(χs)
plot(Ts, χ̄, markersize=10, label="χ")
errorbars!(Ts, χ̄, varχ, whiskerwidth=5)

Cv̄ = mean.(Cvs)
Cv̄ *= maximum(χ̄) / maximum(Cv̄)
varCv = std.(Cvs)
plot!(Ts, Cv̄, markersize=10, label="Cᵥ")
errorbars!(Ts, Cv̄, varCv, whiskerwidth=5)
current_axis().title = "Ising Model behaviour on scale free network with PL exponent 2.0"
current_axis().xlabel = "Temperature T"
axislegend()
save("Ising_CvChi_SF_20.pdf", current_figure())

## phase diagram
function getχ(T; nsamples=100, sz=1000, γ=3.2)
   χs = zeros(nsamples)
   for s in 1:nsamples
      g = scale_free_with_cutoff(sz, γ, 2)

      sv = ones(sz)
      # equilibration for 50 iterations
      simulated_annealing(g, sv; initial_temp=T + 1, end_temp=T, max_steps=100)
      # now with extraction of thermodynamic quantities
      sv, accum = equilibrate(g, sv, T, 100)
      # χ = Var(M)/T
      χs[s] = (accum[4] - accum[3]^2) / T
      # C_v = Var(E)/T^2
      #Cvs[id][s] = (accum[2] - accum[1]^2) / T^2
   end
   return mean(χs)
end  # function getχ

"""
   simulateN(T, nsamples, n, kmin, α)
   simulate N samples the ising model on scale free network with the given parameters and
   determine the phase transition temperature from that using iterative refinement
"""
function simulateN(T, nsamples, n, kmin, α)
   #@show T
   mags = zeros(nsamples)
   lck = ReentrantLock()
   @sync Threads.@threads for s in 1:nsamples
      #println("$s")
      g = scale_free_with_cutoff(n, α, kmin)
      sv = ones(n)
      # equilibration for 50 iterations
      simulated_annealing(g, sv; initial_temp=T + 1, end_temp=T, max_steps=100)
      # now with extraction of thermodynamic quantities
      sv, accum = equilibrate(g, sv, T, 100)
      #@show accum
      #Es[s] = accum[1]
      @lock lck begin
         mags[s] = accum[3]
         # χ = Var(M)/T
         #χs[s] = (accum[4] - accum[3]^2) / T
         # C_v = Var(E)/T^2
         #Cvs[s] = (accum[2] - accum[1]^2) / T^2
      end
   end
   mean(abs.(mags)), std(mags)
end  # function simulateN

"""
    myrangelen(r)
    return span of the range in number space
"""
myrangelen(r) = abs(-(extrema(r)...))



function getTransitionTemp(γ=3.2, kmin=1, n=1000; nsamples=100, threshhold=1.0, preferredmag=0.5)
   r = 0:1:5
   nelem = length(r)
   Tcrit = 0.0
   rlength = myrangelen(r)
   evalfunc = Base.Fix2((x, y) -> (x - y)^2, preferredmag)
   while rlength > threshhold
      mags = map(T -> simulateN(T, nsamples, n, kmin, γ)[1], r)
      plot!(collect(r), mags)
      bestmag, bestTpos = findmin(evalfunc, mags)
      bestT = r[bestTpos]
      lines!([bestT, bestT], [0.0, 1.0])
      if bestT == last(r)
         r = r .+ rlength / 2
      elseif bestT == first(r) && bestT > 0
         r = r .- min(rlength / 2, first(r))
      else
         newhalflength = rlength / (nelem / 2)
         r = range(max(bestT - newhalflength, 0), bestT + newhalflength, length=nelem)
         rlength = myrangelen(r)
         Tcrit = bestT
      end
   end
   Tcrit
end

using CairoMakie
fig = Figure();
ax = Axis(fig[1, 1]);
fig;

#@time getTransitionTemp(3.0, 2, 500; nsamples=20)
n = 500
kmin = 2
nsamples = 50
γs = collect(2.0:0.2:5.0)


kmin = 2
function sampledtheory(γ, n=n, kmin=kmin, nsamples=nsamples * 4)
   gs = [scale_free_with_cutoff(n, γ, kmin) for _ = 1:nsamples]
   degs = degree.(gs)
   gs = 0
   d2s = mean.([d .^ 2 for d in degs])
   degs = mean.(degs)
   βs = log(1 - mean(degs ./ d2s))
   return -1 / βs
end  # function sampledtheory

theoTs = @showprogress map(γ -> sampledtheory(γ), γs)

importantfig, importantax, lineplt = lines(γs, theoTs, label="theoretical prediction")

current_figure!(fig)
for preferredmag in vcat(0.01:0.01:0.04, 0.05:0.05:0.3)
   @info "now at preferred magnetization $preferredmag"
   transitionTs = @showprogress map(γ -> getTransitionTemp(γ, kmin, n;
         nsamples, threshhold=0.2, preferredmag),
      γs)
   plot!(importantax, γs, transitionTs, label="Target M = $preferredmag")
end

save("Magnetization_progression.pdf", fig)


current_figure!(importantfig)
importantfig
yl = "Critical temperature"
xl = "Powerlaw exponent γ"
importantax.title = "$yl vs $xl"
importantax.xlabel = xl
importantax.ylabel = yl
axislegend()
save("phase_diagram_cutoff_variable_powerlaw.pdf", importantfig)

## erdos_renyi
# size of networks
n = 500

# probability of edge formations
ps = [1 / 500, 1 / 250, 1 / 100, log(500) / 500,
   0.015, 0.02, 0.04, 0.06, 0.08, 0.1]

# temperatures to simulate ising model for
Ts = 10 .^ (-1:0.1:2.0)


# number of networks to draw from the ensemble
nemsemble = 100

mags = zeros(nemsemble, length(Ts), length(ps))

lck = ReentrantLock()
@showprogress for p_i in eachindex(ps)
   p = ps[p_i]
   for (id, T) in enumerate(Ts)
      Threads.@threads for s in 1:nemsemble
         g = erdos_renyi(n, p)
         sv = ones(nv(g))
         # equilibration for 100 iterations
         simulated_annealing(g, sv; initial_temp=T + 1, end_temp=T, max_steps=100)
         # now with extraction of thermodynamic quantities
         sv, accum = equilibrate(g, sv, T, 100)
         @lock lck (mags[s, id, p_i] = accum[3])
      end
      GC.gc()
   end
end

fig = Figure(; size=(700, 350))
ax = Axis(fig[1, 1])

meanmag = abs.(mean(mags; dims=1))
varmag = std(mags; dims=1)

colors = CairoMakie.Colors.distinguishable_colors(10, [colorant"white"], dropseed=true)

for p_i in eachindex(ps)
   thismeanmag = meanmag[1, :, p_i]
   p = round(ps[p_i], sigdigits=3)
   kmean = round((n-1)*p, sigdigits=3)
   plot!(Ts, thismeanmag, markersize=10, label="M for p=$p => k̄=$kmean", color=colors[p_i])
   errorbars!(Ts, thismeanmag, varmag[1, :, p_i], whiskerwidth=5, label="M for p=$p => k̄=$kmean", color=colors[p_i])
end

ax.xscale=log
ax.title = "Ising Model behaviour on Erdos-Renyi networks with various degrees"
ax.ylabel = "Magnetization per Vertex M"
ax.xlabel = "Temperature T"
fig[1, 2] = Legend(fig, ax, merge=true)

save("Ising_magnetization_ER_many.pdf", fig)


## watts-strogatz
# size of networks
n = 500

# probability of edge rewiring
ps = collect(0:0.1:1)

# temperatures to simulate ising model for
Ts = collect(0:0.2:5)

# number of initial neighbours
k = 4

# number of networks to draw from the ensemble
nemsemble = 100

mags = zeros(nemsemble, length(Ts), length(ps))

lck = ReentrantLock()
@showprogress for p_i in eachindex(ps)
   p = ps[p_i]
   for (id, T) in enumerate(Ts)
      Threads.@threads for s in 1:nemsemble
         g = watts_strogatz(n, k, p)
         sv = ones(nv(g))
         # equilibration for 100 iterations
         simulated_annealing(g, sv; initial_temp=T + 1, end_temp=T, max_steps=100)
         # now with extraction of thermodynamic quantities
         sv, accum = equilibrate(g, sv, T, 100)
         @lock lck (mags[s, id, p_i] = accum[3])
      end
      GC.gc()
   end
end

fig = Figure(; size=(700, 350))
ax = Axis(fig[1, 1])

meanmag = abs.(mean(mags; dims=1))
varmag = std(mags; dims=1)

colors = CairoMakie.Colors.distinguishable_colors(length(ps), [colorant"white"], dropseed=true)

for p_i in eachindex(ps)
   thismeanmag = meanmag[1, :, p_i]
   p = round(ps[p_i], sigdigits=3)
   kmean = round((n-1)*p, sigdigits=3)
   plot!(Ts, thismeanmag, markersize=10, label="M for p=$p", color=colors[p_i])
   errorbars!(Ts, thismeanmag, varmag[1, :, p_i], whiskerwidth=5, label="M for p=$p", color=colors[p_i])
end

ax.title = "Ising Model behaviour on Watts-Strogatz networks with various rewiring probabilities"
ax.ylabel = "Magnetization per Vertex M"
ax.xlabel = "Temperature T"
fig[1, 2] = Legend(fig, ax, merge=true)

save("Ising_magnetization_WS_many.pdf", fig)


## SBM
# size of networks
ncomms = 8
percomm = 125
n = fill(percomm, ncomms)

# probability of intra community edges
intrap = 0.05

# relative probability of intercommunity edges
ps = [0, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 0.2, 0.5, 0.8]

# temperatures to simulate ising model for
Ts = collect(range(0, 45, 50))

# number of networks to draw from the ensemble
nemsemble = 100

mags = zeros(nemsemble, length(Ts), length(ps))

lck = ReentrantLock()
@showprogress for p_i in eachindex(ps)
   p = ps[p_i]
   for (id, T) in enumerate(Ts)
      Threads.@threads for s in 1:nemsemble
         sbm_matrix = probsTomeanDeg(intrap, intrap*p, n)
         g = stochastic_block_model(sbm_matrix, n)
         sv = ones(nv(g))
         # equilibration for 100 iterations
         simulated_annealing(g, sv; initial_temp=T + 1, end_temp=T, max_steps=100)
         # now with extraction of thermodynamic quantities
         sv, accum = equilibrate(g, sv, T, 100)
         @lock lck (mags[s, id, p_i] = accum[3])
      end
      GC.gc()
   end
end

fig = Figure(; size=(700, 350))
ax = Axis(fig[1, 1])

meanmag = abs.(mean(mags; dims=1))
varmag = std(mags; dims=1)

colors = CairoMakie.Colors.distinguishable_colors(length(ps), [colorant"white"], dropseed=true)

for p_i in eachindex(ps)
   thismeanmag = meanmag[1, :, p_i]
   p = round(ps[p_i], sigdigits=3)
   plot!(Ts, thismeanmag, markersize=10, label="M for p_ext/p_int=$p", color=colors[p_i])
   errorbars!(Ts, thismeanmag, varmag[1, :, p_i], whiskerwidth=5, label="M for p_ext/p_int=$p", color=colors[p_i])
end

ax.title = "Ising Model behaviour on Stochastic block model networks\nwith various relative inter community connections probabilities"
ax.ylabel = "Magnetization per Vertex M"
ax.xlabel = "Temperature T"
fig[1, 2] = Legend(fig, ax, merge=true)

save("Ising_magnetization_SBM_many.pdf", fig)
