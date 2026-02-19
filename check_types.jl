
include("CSRGraph.jl")
include("Dijkstra.jl")
include("DMMSY-SSSP.jl")

using .CSRGraphModule
using .DMMSYSSSP
using InteractiveUtils

n, m = 1000, 5000
g = random_graph(n, m, 100.0)
println("--- Type Stability Check: ssp_duan ---")
@code_warntype ssp_duan(g, 1)
