#!/usr/bin/env julia

include("CSRGraph.jl")
include("Dijkstra.jl")
include("DMMSY-SSSP.jl")
include("DMMSY_Research.jl")

using .CSRGraphModule
using .DijkstraModule
using .DMMSYSSSP
using .DMMSYResearch
using Printf
using Random

function run_full_benchmark()
    println("SSSP Unified Benchmark")
    println("="^74)
    
    test_cases = [
        (1000, 5000, 100),
        (10000, 50000, 50),
        (25000, 125000, 30),
        (50000, 250000, 20),
        (100000, 500000, 15),
        (125000, 625000, 12),
        (150000, 750000, 10),
        (175000, 875000, 8),
        (200000, 1000000, 8),
        (300000, 1500000, 5),
        (400000, 2000000, 4),
        (500000, 2500000, 4),
        (750000, 3750000, 2),
        (1000000, 5000000, 2)
    ]

    @printf("%8s %8s %12s %12s %12s %8s %8s\n", "n", "m", "Dijkstra", "DMMSY Opt", "DMMSY Res", "SpdOpt", "SpdRes")
    println("-"^74)

    for (n, m, trials) in test_cases
        g = random_graph(n, m, 100.0)
        source = 1

        # Warm-up
        DijkstraModule.dijkstra_ref(g, source)
        DMMSYSSSP.ssp_duan(g, source)
        DMMSYResearch.ssp_duan_research(g, source)

        t_dij = @elapsed for _ in 1:trials DijkstraModule.dijkstra_ref(g, source) end
        t_opt = @elapsed for _ in 1:trials DMMSYSSSP.ssp_duan(g, source) end
        t_res = @elapsed for _ in 1:trials DMMSYResearch.ssp_duan_research(g, source) end

        t_dij /= trials
        t_opt /= trials
        t_res /= trials

        spd_opt = t_dij / t_opt
        spd_res = t_dij / t_res

        @printf("%8d %8d %12.4fms %12.4fms %12.4fms %8.2fx %8.2fx\n", 
                n, m, t_dij*1000, t_opt*1000, t_res*1000, spd_opt, spd_res)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_full_benchmark()
end
