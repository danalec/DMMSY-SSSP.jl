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
        (100, 500, 20),
        (500, 2000, 20),
        (1000, 5000, 15),
        (2500, 12500, 10),
        (5000, 25000, 5),
        (10000, 50000, 3),
        (25000, 125000, 2),
        (50000, 250000, 1)
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

            # Warm-up (Ensures JIT compilation is finished before recording)
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
