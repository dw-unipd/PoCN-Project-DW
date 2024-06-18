## imports
using Graphs, MetaGraphsNext
using Glob
using CSV, DataFrames
using GeometryBasics
using ProgressMeter
using CairoMakie, PairPlots

## files
ispath("code/") && cd("code/task_41/")
nodefiles = glob("*.nodes", "../../data")
edgefiles = glob("*.edges", "../../data")
@show length(nodefiles) == length(edgefiles)

cities = [splitext(f)[1] for f in nodefiles]
ncities = length(cities)

## data types, adjusted for the analysis
# only care about where and how called
struct nodeData
   name::AbstractString
   lonlat::Point2
end

struct edgeData
   length::Float64
end

## helper functions
include("common_functions.jl")
include("../common_utils.jl")

# custom distance function for use in the graph constructor
function graphDistance(edge)
   # edge is either nothing for interline connections or an instance of edgeData
   isnothing(edge) && return 0.0
   edge.length
end  # function graphDistance

common_parse_args = (delim=';', comment="#")

function readGraph(city)
   g = MetaGraph(Graph(), Int, nodeData, Union{edgeData,Nothing},
      "Graph for $(last(splitpath(city)))", graphDistance, Inf)
   node_data = CSV.read(city * ".nodes", DataFrame; common_parse_args..., types=[Int, String, Float64, Float64, Int, Int])
   edge_data = CSV.read(city * ".edges", DataFrame; common_parse_args..., types=[Int, Int, Int, String, Int])
   sort!(node_data, :id)
   for row in eachrow(node_data)
      name = row.name
      if ismissing(name)
         name = ""
      end
      if ismissing(row.id)
         @warn "missing $row"
         @error city
      end
      g[row.id] = nodeData(name, Point2([row.longitude, row.latitude]))
      #println(g[row.id])
   end
   for row in eachrow(edge_data)
      from = row.from_id
      to = row.to_id
      distance = try
         dist(
            g[from].lonlat, g[to].lonlat
         )
      catch e
         @warn e
         @info city
         0.0
      end
      #println("$city: $from -> $to = $distance")
      g[from, to] = edgeData(distance)
   end
   g
end  # function readGraph

"""
   connectGraph(g, id2name)
   g: Graph
   id2name: mapping station id to its name, faster lookup than using the graph
   returns graph with duplicate stations removed
   if the first part of the station name of node N matches an already encountered one N0,
      all the neighbours of N are connected to N0 and N is removed from the graph
"""
function connectGraph(g, id2name; minimum_match=20, connect_dist=100)
   name2id = Dict()
   while true
      oldl = length(labels(g))
      for label in labels(g)
         fullname = id2name[label]
         # no sense in trying to match multiple stations without name
         isempty(fullname) && continue
         locname = ""
         for str in split(fullname, ' ')
            locname *= (isempty(locname) ? "" : ' ') * str
            length(locname) > minimum_match && break
         end
         repl_id = get!(name2id, locname, label)
         label == repl_id && continue
         locname != fullname && println("removing duplicate $locname ($fullname)")
         neighbours = neighbor_labels(g, label)
         for neighbour in neighbours
            g[repl_id, neighbour] = g[label, neighbour]
         end
         rem_vertex!(g, code_for(g, label))
         break
      end
      length(labels(g)) == oldl && break
   end
   # connect nearby stations
   for l1 in labels(g)
      for l2 in labels(g)
         l1 == l2 && continue
         haskey(g, l1, l2) && continue
         distance = dist(
            g[l1].lonlat, g[l2].lonlat
         )
         if distance < connect_dist
            g[l1, l2] = edgeData(distance)
         end
      end
   end
   # remove self loops
   for label in labels(g)
      code = code_for(g, label)
      rem_edge!(g, code, code)
   end
   g
end  # function canonicalizeGraph

## data per city
# collection of properties of interest to compare across cities

unconnected_df = DataFrame(;
   # network size
   sizes=zeros(ncities),
   # mean degree
   meandegs=zeros(ncities),
   # degree assortativity
   assortativities=zeros(ncities),
   # diameter in meters
   diameters_real=zeros(ncities),
   # diameter in # of edges
   diameters_edges=zeros(ncities),
   # size lcc
   lccs=zeros(ncities),
   # powerlaw exponent
   αs=Union{Missing,Float64}[zeros(ncities)...]
)

connected_df = DataFrame(;
   sizes=zeros(ncities),
   meandegs=zeros(ncities),
   assortativities=zeros(ncities),
   diameters_real=zeros(ncities),
   diameters_edges=zeros(ncities),
   lccs=zeros(ncities),
   αs=Union{Missing,Float64}[zeros(ncities)...]
)

## main processing loop

@showprogress for i in eachindex(cities)
   city = cities[i]
   cityG = readGraph(city)
   unconnected_df[i, :sizes] = nv(cityG)
   unconnected_df[i, :meandegs] = mean(degree(cityG))
   unconnected_df[i, :assortativities] = begin
      ass = assortativity(cityG)
      isnan(ass) ? 0.0 : ass
   end
   lcc = argmax(length, connected_components(cityG))
   unconnected_df[i, :lccs] = length(lcc) / nv(cityG)
   subG, indexmap = induced_subgraph(cityG, lcc)
   unconnected_df[i, :diameters_real] = diameter(subG)
   unconnected_df[i, :diameters_edges] = diameter(subG, Graphs.DefaultDistance())
   # the powerlaw fit can fail if igraph detects a negative alpha, this is captured here
   powerlaw_result = try
      getPowerLawFit(cityG)
   catch e
      if e isa PyCall.PyError
         missing
      else
         rethrow(e)
      end
   end
   unconnected_df[i, :αs] = if ismissing(powerlaw_result) || powerlaw_result.p < 0.05
      missing
   else
      powerlaw_result.alpha
   end
   connG = connectGraph(deepcopy(cityG), Dict((l => cityG[l].name for l in labels(cityG))))
   connected_df[i, :meandegs] = mean(degree(connG))
   connected_df[i, :sizes] = nv(connG)
   connected_df[i, :assortativities] = begin
      ass = assortativity(connG)
      isnan(ass) ? 0.0 : ass
   end
   lcc = argmax(length, connected_components(connG))
   connected_df[i, :lccs] = length(lcc) / nv(connG)
   subG, indexmap = induced_subgraph(connG, lcc)
   connected_df[i, :diameters_real] = diameter(subG)
   connected_df[i, :diameters_edges] = diameter(subG, Graphs.DefaultDistance())
   # the powerlaw fit can fail if igraph detects a negative alpha, this is captured here
   powerlaw_result = try
      getPowerLawFit(connG)
   catch e
      if e isa PyCall.PyError
         missing
      else
         rethrow(e)
      end
   end
   connected_df[i, :αs] = if ismissing(powerlaw_result) || powerlaw_result.p < 0.05
      missing
   else
      powerlaw_result.alpha
   end
end

## plot both 

fig = Figure(; size=(1200, 1200))
pairplot(fig[1, 1], unconnected_df => (
   PairPlots.Scatter(; markersize=5),
   PairPlots.MarginHist(),
   PairPlots.MarginDensity(),
   PairPlots.MarginConfidenceLimits()
))
Label(fig[0, :],
   "Statistics of networks as extracted from the data",
   tellwidth=false)
fig
save("pairplot_unconnected.png", fig)

fig = Figure(; size=(1200, 1200))
pairplot(fig[1, 1], connected_df => (
   PairPlots.Scatter(; markersize=5),
   PairPlots.MarginHist(),
   PairPlots.MarginDensity(),
   PairPlots.MarginConfidenceLimits(),
))
Label(fig[0, :],
   "Statistics of networks after adding probable connections",
   tellwidth=false)
current_figure()
save("pairplot_connected.png", fig)

## plot both
fig = Figure(; size=(1200, 1200))
pairplot(fig[1, 1],
   PairPlots.Series(unconnected_df,
      label="as extracted",
      color=Makie.wong_colors(0.6)[3]) => (
      PairPlots.Scatter(markersize=4),
      #PairPlots.Contourf(),
      PairPlots.MarginDensity(
         linewidth=2.5f0
      ),
   ),
   PairPlots.Series(connected_df,
      label="after connecting",
      color=Makie.wong_colors(0.6)[2]) => (
      PairPlots.Scatter(markersize=4),
      #PairPlots.Contour(),
      PairPlots.MarginDensity(
         linewidth=2.5f0
      ),
   ),
   axis=(;
      sizes=(; scale=log10),
      meandegs=(scale=log10,),
      diameters_real=(scale=log10,),
      diameters_edges=(scale=log10,),
   ),
   diagaxis=(;
      xgridvisible=true,
   ),
   bodyaxis=(;
      ygridvisible=true,
      xgridvisible=true,
   )
)
Label(fig[0, :],
   "Statistics of networks compared between graph families",
   tellwidth=false)
save("../../latex/images/pairplot_both.pdf", fig)
fig

fig, ax, plt = plot(unconnected_df.αs, unconnected_df.assortativities)
plot!(connected_df.αs, connected_df.assortativities)
ax.xscale=log10
fig