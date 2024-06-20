using Graphs
using PyCall
ig = pyimport("igraph")
function ig2jg(ig_g)
   jg_g = SimpleGraph(ig_g.vcount())
   for e in ig_g.es()
       add_edge!(jg_g, e.source + 1, e.target + 1)
   end
   return jg_g
end

function jg2ig(jg_g)
   nvert = nv(jg_g)
   es = [Tuple(e).-1 for e in edges(jg_g)]
   return ig.Graph(nvert, es)
end  # function jg2ig

function getPowerLawFit(jg)
   igg = jg2ig(jg)
   data = [x[3] for x in igg.degree_distribution().bins()]
   ig.power_law_fit(data)
end  # function getPowerLawFit