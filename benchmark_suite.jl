#!/usr/bin/env julia
"""
benchmark_suite.jl — Performance Data Gatherer for DMMSY-SSSP

This script runs a comprehensive benchmark across various graph sizes and
saves the weighted average timing results to 'benchmark_data.csv'.
"""

include("CSRGraph.jl")
include("Dijkstra.jl")
include("DMMSY-SSSP.jl")
include("DMMSY_Research.jl")

using .CSRGraphModule
using .DijkstraModule
using .DMMSYSSSP
using .DMMSYResearch
using Printf

function run_benchmark_collection()
    # Configuration: (n, m, trials)
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

    println("DMMSY Data Collection Started...")
    println("="^40)

    output_file = "benchmark_data.csv"
    
    open(output_file, "w") do f
        # Write CSV Header
        println(f, "n,m,dijkstra,dmmsy_opt,dmmsy_res")
        
        for (n, m, trials) in test_cases
            print("  Processing n=$n, m=$m ($trials trials)... ")
            
            g = random_graph(n, m, 100.0)
            source = 1

            # Robust Warm-up
            g_warm = random_graph(100, 500, 100.0)
            for _ in 1:5
                DijkstraModule.dijkstra_ref(g_warm, 1)
                DMMSYSSSP.ssp_duan(g_warm, 1)
                DMMSYResearch.ssp_duan_research(g_warm, 1)
            end
            
            # Specific Warm-up
            DijkstraModule.dijkstra_ref(g, source)
            DMMSYSSSP.ssp_duan(g, source)
            DMMSYResearch.ssp_duan_research(g, source)

            # Execution
            t_dij = @elapsed for _ in 1:trials DijkstraModule.dijkstra_ref(g, source) end
            t_opt = @elapsed for _ in 1:trials DMMSYSSSP.ssp_duan(g, source) end
            t_res = @elapsed for _ in 1:trials DMMSYResearch.ssp_duan_research(g, source) end

            # Averaging (convert to ms)
            avg_dij = (t_dij / trials) * 1000
            avg_opt = (t_opt / trials) * 1000
            avg_res = (t_res / trials) * 1000

            # Log to file
            @printf(f, "%d,%d,%.6f,%.6f,%.6f\n", n, m, avg_dij, avg_opt, avg_res)
            
            @printf("dijkstra: %.4fms, opt: %.4fms, res: %.4fms -> Spd: %.2fx\n", avg_dij, avg_opt, avg_res, avg_dij/avg_opt)
            println("Done.")
        end
    end

    println("="^40)
    println("Data successfully saved to: $output_file")
    println("You can now copy the CSV values into 'benchmark_results.html'.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark_collection()
end
