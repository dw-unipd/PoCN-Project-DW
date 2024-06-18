struct nodeData
   name::AbstractString
   lonlat::Point2
   mode::Int
   year::Int
end

struct edgeData
   mode::Int
   line::AbstractString
   year::Int
   length::Float64
end