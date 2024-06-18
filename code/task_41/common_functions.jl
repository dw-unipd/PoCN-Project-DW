using Distances

# custom distance function for use in the graph constructor
function graphDistance(edge)
   # edge is either nothing for interline connections or an instance of edgeData
   isnothing(edge) && return 0.0
   edge.length
end  # function graphDistance

h = Haversine()
dist(x, y) = h(x, y) # geodesic distance between 2 points
dist(tup::Tuple) = dist(tup...)
dist(x::Union{Point,Real}) = y -> dist(x, y) # create function to use for map over collection
