## imports
# reading data
using CSV
# general data structure
using DataFrames
# for network creation and analysis
using Graphs, MetaGraphsNext
# for mathematical vectors, convenient distance calculation
using GeometryBasics, Distances, LinearAlgebra
# mainly to find stations near tracks
using NearestNeighbors
# plotting libraries
using GraphMakie, CairoMakie
# functions such as the mean come from here
using Statistics
# for status report on running computation
using ProgressMeter

## helper structs
include("common_types.jl")

## utility functions
include("common_functions.jl")

# position strings are wrapped in TYPEDECLARATION()
# this function extracts the numbers between the parentheses
function extractGeometry(str)
   try
      startpos = findfirst('(', str)
      endpos = findlast(')', str)
      any(isnothing, [startpos, endpos]) && return nothing
      str[startpos+1:endpos-1]
   catch e
      @show str
      rethrow(e)
   end
end # function extractGeometry

# convert POINT(...) to Point2 datastructure
function getPoint(str)
   coordstr = extractGeometry(str)
   isnothing(coordstr) && return nothing
   point(coordstr)
end  # function getPoint

# convert LINESEGMENT(...) to vector of Point2
function getPoints(linestr)
   coordstrs = extractGeometry(linestr)
   isnothing(coordstrs) && return nothing
   point.(split(coordstrs, ','))
end  # function getPoints

# each point is "lon lat"
function point(coordstr)
   Point(parse.(Float64, split(coordstr, " "))...)
end  # function point  

function getDistance(s, e, p)
   # distance of Point p to line segment between start s and end e
   l2 = dot(e - s, e - s)
   l2 == 0 && return dist(s, p) # start == end
   # segment parametrized by s + t * (e-s) for t ∈ [0, 1]
   # project s->p onto s->e, clamp to valid t range
   t = clamp(dot(p - s, e - s) / l2,
      0.0, 1.0)
   proj = s + t * (e - s)
   dist(p, proj)
end  # function getDistance

function getDistances(s, e, ps)
   map(p -> getDistance(s, e, p), ps)
end  # function getDistances


# function sort a list of lists of points to their most probable order
function sortSegments(segments)
   while true
      head2heads = map(x -> any(y -> 0 < dist(y, first(x)) < dist(y, last(x)), first.(segments)), segments)#=.|
      map(x -> any(y -> 0 < dist(y, last(x)) < dist(y, first(x)), last.(segments)), segments)=#
      toflip = findall(head2heads)
      if isnothing(toflip) || isempty(toflip)
         break
      else
         nowflip = rand(toflip)
         reverse!(segments[nowflip])
      end
   end

   startseg = findfirst(!,
      # for each segment, check if there is any other segment whose first element is "our" last
      map(x -> any(==(last(x)), first.(segments)),
         segments)
   )


   if isnothing(startseg)
      startseg = findfirst(!,
         # for each segment, check if there is any other segment whose last element is "our" last
         map(x -> any(==(last(x)), last.(segments)),
            segments)
      )
      # if there is, reverse that segment
      isnothing(startseg) || reverse!(segments[startseg])
   end

   isnothing(startseg) ||
      ((segments[startseg], segments[1]) = (segments[1], segments[startseg]))

   # this loop finds the segment that is the most probable continuation
   # --> should come at index idx
   for idx in eachindex(segments)[2:end]
      prevlast = last(segments[idx-1])
      dists = map(dist(prevlast),
         first.(segments[idx:end]))
      #println(round.(dists, sigdigits=2))
      mdist, continuation = findmin(dists)
      do_reverse = false
      if !(mdist ≈ 0.0)
         revdists = map(dist(prevlast),
            last.(segments[idx:end]))
         revmdist, rev_continuation = findmin(revdists)
         if revmdist < mdist
            continuation = rev_continuation
            do_reverse = true
         end
      end
      continuation += idx - 1
      do_reverse && reverse!(segments[continuation])
      continuation == idx ||
         ((segments[idx], segments[continuation]) = (segments[continuation], segments[idx]))
   end
   return segments
end

function sortSegments2(segments)
   sortedseg = first(segments)
   deleteat!(segments, 1)
   while !isempty(segments)
      # start of sorted sequence
      firstfirst = first(sortedseg)
      # minimum distance of start to start of other segments
      fffdist, fffidx = findmin(dist(firstfirst), first.(segments))
      # minimum distance of start to end of other  
      ffldist, fflidx = findmin(dist(firstfirst), last.(segments))
      # find lowest distance and its index and whether it needs to be put in reversed
      ffdist = min(ffldist, fffdist)
      ffidx = ffldist < fffdist ? fflidx : fffidx
      # should reverse if start-end distance is smaller than start-start
      # >--> <--< should be reversed
      ffreverse = ffldist > fffdist

      # do same things for end of sequence
      # end of sorted segments
      lastlast = last(sortedseg)
      llfdist, llfidx = findmin(dist(lastlast), first.(segments))
      llldist, lllidx = findmin(dist(lastlast), last.(segments))
      lldist = min(llfdist, llldist)
      llidx = llfdist < llldist ? llfidx : lllidx
      llreverse = llldist < llfdist

      # combine knowledge
      # insert at front: true, append at end: false
      shouldinsert = ffdist < lldist
      shouldreverse = shouldinsert ? ffreverse : llreverse
      minidx = shouldinsert ? ffidx : llidx
      toadd = segments[minidx]
      shouldreverse && reverse!(toadd)
      if shouldinsert
         sortedseg = vcat(toadd, sortedseg)
      else
         sortedseg = vcat(sortedseg, toadd)
      end
      deleteat!(segments, minidx)
   end # while
   sortedseg
end  # function sortSegments2

function remove_duplicates(series)
   isempty(series) && return series
   out = eltype(series)[first(series)]
   sizehint!(out, length(series))
   prevval = first(out)
   for value in @view series[2:end]
      value == prevval && continue
      push!(out, value)
      prevval = value
   end
   out
end  # function remove_duplicates

## data loading
ispath("code/") && cd("code/task_41/")
# load the data, data comes from citylines.com
# each file also contains data related to the web interface, those are not further described

# details of available cities
# id, name, center coordinates, start of Transport 
# country and sometimes state, length of network
cities = CSV.read("data/cities.csv", DataFrame)

# details of available lines
# id, city it belongs to, name, system it belongs, transport mode
lines = CSV.read("data/lines.csv", DataFrame)

# details of available stations
# id, name, geographical location, buildstart, opening an closure year, city id
stations = CSV.read("data/stations.csv", DataFrame)

# mapping between stations and lines, info on when this part was/is in use
stationlines = CSV.read("data/stationlines.csv", DataFrame)

# mapping between mode ids and their descriptions
modes = CSV.read("data/modes.csv", DataFrame)

# details on mapped sections, segments or complete lines of transport
# id, list of waypoints, time information, length of section, city id
sections = CSV.read("data/sections.csv", DataFrame)

# mapping between sections and line ids
# operational time frame also listed
sectionlines = CSV.read("data/sectionlines.csv", DataFrame)

# details on available systems, e.g. metro and tram are different systems in one city
# id, city id, name, whether current or historic system, length of mapped tracks
systems = CSV.read("data/systems.csv", DataFrame)


## main processing area
p = Progress(nrow(systems); showspeed=true)

# assume stations are never more than this far from their respective line
dist_thresh = 20 #meters

@time for city in eachrow(cities)
   # for each city extract the relevant parts of the total DataFrames
   # to speed up computations
   # the city_id is the key
   city_id = city.id
   city_name = city.name
   city_lines = lines[lines.city_id.==city_id, :]
   city_stations = stations[stations.city_id.==city_id, :]
   city_stationlines = stationlines[stationlines.city_id.==city_id, :]
   city_sections = sections[sections.city_id.==city_id, :]
   city_sectionlines = sectionlines[sectionlines.city_id.==city_id, :]
   city_systems = systems[systems.city_id.==city_id, :]

   # can't build a graph for an empty city
   if (nrow(city_stations) == 0
       || nrow(city_sections) == 0
       || nrow(city_lines) == 0)
      continue
   end

   default_year = city.start_year

   g = MetaGraph(Graph(), Int, nodeData, Union{edgeData,Nothing}, "Graph for $(city.name)", graphDistance)
   fig = Figure()
   ax = Axis(fig[1, 1])

   # treat each system separately
   for system in eachrow(city_systems)
      # update progress meter
      next!(p)
      system_id = system.id

      system_lines = city_lines[city_lines.system_id.==system_id, :]

      # process each line separately
      for line = eachrow(system_lines)
         line_id = line.id
         line_stationlines = city_stationlines[city_stationlines.line_id.==line_id, :]
         line_stations = city_stations[in.(city_stations.id, (line_stationlines.station_id,)), :]

         # if no stations are here, nothing can contribute to the graph
         iszero(nrow(line_stations)) && continue

         line_sectionlines = city_sectionlines[city_sectionlines.line_id.==line_id, :]
         line_sections = city_sections[in.(city_sections.id, (line_sectionlines.section_id,)), :]


         # if no sections are here, nothing can contribute to the graph
         iszero(nrow(line_sections)) && continue

         # get coordinates of all the stations
         station_coords = getPoint.(line_stations.geometry)

         # all connections should be shorter than this
         # as distances tend to both get longer and spread out with distance from center
         maxdist = begin
            coords = getPoint.(unique(line_stations, :name).geometry)
            # get minimum distance of all stations to any other of a different name
            # provides a reasonable upper bound on the length of connections
            mindists = map(x -> minimum(y -> y == x ? Inf : dist(x, y), coords), coords)
            5 * maximum(mindists)
         end
         station_ids = line_stations.id
         # create tree to find nearest neighbours of 
         nntree = BallTree(station_coords, Haversine())

         # list of segments
         # each segment is a list of points that signify part of that line
         line_segments = getPoints.(line_sections.geometry)
         filter!(!isnothing, line_segments)

         # if no segments, nothing can be won here
         isnothing(line_segments) && continue
         isempty(line_segments) && continue

         # the sort implementation can sometimes take multiple runs to converge
         #@enter 
         sorted_segments = sortSegments2(deepcopy(line_segments))
         #=while sorted_segments != line_segments
            line_segments = sorted_segments
            sorted_segments = sortSegments(deepcopy(line_segments))
         end=#

         # join all segments into one list, remove coordinates occuring more than once in a row
         line_points = remove_duplicates(sorted_segments)

         plot!(ax, line_points, markersize=5)
         plot!(ax, station_coords, marker=:x)
         # get the closest stations to each line_point
         ids, dists = knn(nntree, line_points, min(5, length(nntree.data)), true)

         # list of stations in order of connection, contains the id into the station_ids/station_coords/line_stations
         # so station_ids can be indexed by the entries 
         connected_stations = Int[]

         dists[1][1] < dist_thresh && push!(connected_stations, ids[1][1])

         for (idx, (s, e)) in enumerate(zip(line_points[1:end-1], line_points[2:end]))
            local stations = ids[idx]
            line_dist, stationid = findmin(station->getDistance(s, e, station_coords[station]), stations)
            station = stations[stationid]
            # don't want to duplicate station
            isempty(connected_stations) || station == last(connected_stations) && continue
            # station is too far from line at this point
            line_dist > dist_thresh && continue
            push!(connected_stations, station)
         end

         # map all the ids to the first station with the same name
         # remove duplicate stations, e.g. different directions
         map!(connected_stations, connected_stations) do idx
            findfirst(==(line_stations.name[idx]), line_stations.name)
         end
         connected_stations = remove_duplicates(connected_stations)

         # create graph vertices for all connected stations
         for idx in unique(connected_stations)
            row = line_stations[idx, :]
            year = row.opening
            # if the year for this data point is missing, use the one from the city as a fallback
            year = (ismissing(year) || iszero(year) || isnothing(year)) ? default_year : year
            # remove newlines from name
            name = replace(strip(row.name), "\n"=>" ")
            g[row.id] = nodeData(name, station_coords[idx], line.transport_mode_id, year)
         end

         # add edges for successive entries in connected_stations
         for (idx, dest) in pairs(station_ids[connected_stations[2:end]])
            # idx is shifted by one vs entry of connected_stations -> perfectly what we want
            source = station_ids[connected_stations[idx]]
            # have already added that edge in aprevious iteration
            haskey(g, source, dest) && continue
            edge_length = dist(g[source].lonlat, g[dest].lonlat)
            # not a reasonable connection if longer than above calculated threshhold
            edge_length > maxdist && continue
            # assume all edges are not younger than their stations
            # possibly not perfect but saves *a lot* of computation and should be a reasonable assumption
            year = max(g[source].year, g[dest].year)
            g[source, dest] = edgeData(line.transport_mode_id, line.name, year, edge_length)
         end
      end # lines
   end # systems

   # this apparently has stations and sections but not together at one place
   # e.g. Berlin
   # analysing a system with only one station or one edge is nonsensical
   if nv(g) < 2 || ne(g) < 2
      continue
   end

   graphplot!(ax, g, layout=[g[l].lonlat for l in labels(g)],
      nodemarker=:+, node_size=7)
   ax.title = "$city_name"
   save("output/$city_name.pdf", fig)

   # dictionary to map internal ids to continous ids starting from 1
   labelmap = Dict(zip(labels(g), 1:nv(g)))

   # write node file
   open("output/$city_name.nodes", "w") do f
      println(f, "# node data for $city_name, $(city.country)")
      println(f, "id;name;latitude;longitude;mode;year")
      for l in labels(g)
         data = Any[labelmap[l]]
         gdata = g[l]
         lon, lat = gdata.lonlat
         append!(data, [gdata.name, lat, lon, gdata.mode, gdata.year])
         println(f, join(data, ';'))
      end
   end

   open("output/$city_name.edges", "w") do f
      println(f, "# edge data for $city_name, $(city.country)")
      println(f, "from_id;to_id;mode;line;year")
      for e in edge_labels(g)
         from, to = e
         data = Any[labelmap[from], labelmap[to]]
         edata = g[from, to]
         # one of the fake edges introduced for connectedness
         isnothing(edata) && continue
         append!(data, [edata.mode, edata.line, edata.year])
         println(f, join(data, ';'))
      end
   end
   #println("done with city $city_name")
end # cities
