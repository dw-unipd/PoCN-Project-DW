using Distributions
include("py_igraph_Graphs_interop.jl")
"""
   scale_free_with_cutoff(n, α, cutoff)
   n: number of nodes
   α: power law exponent
   cutoff: minimum degree of sequence
   inspired by https://stackoverflow.com/questions/15924767/generate-scale-free-network-with-fixed-minimum-degree-using-igraph
"""
function scale_free_with_cutoff(n, α, cutoff)
   dist = Pareto(α, cutoff)
   degrees = round.(Int, rand(dist, n))
   ntries = 1
   while !isgraphical(degrees)
      @debug "going for try $ntries in scale_free_with_cutoff"
      degrees = round.(Int, rand(dist, n))
      ntries += 1
   end
   #random_configuration_model(n, degrees)
   if Threads.threadid() != 1
      if Threads.threadpoolsize(:interactive) == 1
         fetch(Threads.@spawn :interactive ig2jg(ig.Graph.Degree_Sequence(degrees, method="vl")))
      else
         random_configuration_model(n, degrees)
      end
   else
      ig2jg(ig.Graph.Degree_Sequence(degrees, method="vl"))
   end
end


"""
    size_lcc(g)
    g: Graph
    returns size of largest connected component of g
"""
function size_lcc(g)
   nv(g) == 0 && return 0
   maximum(length.(connected_components(g)))
end  # function size_lcc

"""
   probsTomeanDeg(ps::AbstractMatrix, n::AbstractVector)
   igraph takes probabilities, where as Graphs.jl takes average degree of nodes
   this function converts the probability matrix used by igraph::sample_sbm
   to the matrix used by Graphs.stochastic_block_model
"""
function probsTomeanDeg(ps::AbstractMatrix, n::AbstractVector)
   nmat = [(n[max(i, j)] - (i == j)) for i = 1:length(n), j = 1:length(n)]
   [nmat[i] * ps[i] for i in CartesianIndices(ps)]
end  # function probsTomeanDeg



"""
   probsTomeanDeg(p_int<:Real, p_ext<:Real, n::AbstractVector)
   convert a single internal probability and external probability
   to the format taken by Graph.jl's constructor
"""
function probsTomeanDeg(p_int :: Real, p_ext :: Real, n::AbstractVector)
   p_ints = fill(p_int, length(n))
   probsTomeanDeg(Graphs.SimpleGraphs.sbmaffinity(p_ints, p_ext, n), n)
end  # function probsTomeanDeg
