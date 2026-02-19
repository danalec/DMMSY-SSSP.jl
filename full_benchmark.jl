#!/usr/bin/env julia
"""
full_benchmark.jl — High-Precision Performance Data Gatherer
"""

include("CSRGraph.jl")
include("Dijkstra.jl")
include("DMMSY-SSSP.jl")

using .CSRGraphModule
using .DijkstraModule
using .DMMSYSSSP
using Printf

function run_full_benchmark()
    # Configuration: (n, m, trials)
    # Increased trials to reduce noise and show true algorithmic advantage
    # Comprehensive sequence from 1k to 1M
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

    println("DMMSY High-Precision Benchmark Started...")
    
    if !DMMSYSSSP.verify_correctness()
        println("Correctness check FAILED. Aborting benchmark.")
        return
    end
    
    println("="^60)
    @printf("%-10s %-10s %-12s %-12s %-10s\n", "n", "m", "Dijkstra", "DMMSY Opt", "Speedup")
    println("-"^60)

    output_file = "benchmark_data.csv"
    
    open(output_file, "w") do f
        println(f, "n,m,dijkstra,dmmsy_opt,dmmsy_res")
        
        for (n, m, trials) in test_cases
            g = random_graph(n, m, 100.0)
            source = 1

            # Robust Warm-up
            g_warm = random_graph(100, 500, 100.0)
            for _ in 1:10
                DijkstraModule.dijkstra_ref(g_warm, 1)
                DMMSYSSSP.ssp_duan(g_warm, 1)
            end
            
            # Specific Warm-up
            DijkstraModule.dijkstra_ref(g, source)
            DMMSYSSSP.ssp_duan(g, source)

            # Execution
            t_dij = @elapsed for _ in 1:trials DijkstraModule.dijkstra_ref(g, source) end
            t_opt = @elapsed for _ in 1:trials DMMSYSSSP.ssp_duan(g, source) end

            # Averaging (ms)
            avg_dij = (t_dij / trials) * 1000
            avg_opt = (t_opt / trials) * 1000
            spd = avg_dij / avg_opt

            @printf("%-10d %-10d %-10.4fms %-10.4fms %-8.2fx\n", n, m, avg_dij, avg_opt, spd)
            @printf(f, "%d,%d,%.6f,%.6f,%.6f\n", n, m, avg_dij, avg_opt, 0.0)
        end
    end

    println("="^60)
    println("Full benchmark complete. Results saved to $output_file")
    
    # Auto-embed into HTML for offline viewing (bypasses Fetch CORS)
    try
        csv_data = read(output_file, String)
        html_path = "benchmark_results.html"
        if isfile(html_path)
            html_content = read(html_path, String)
            # Find the rawData block and replace it
            new_data_js = "const EMBEDDED_DATA = `$(csv_data)`;"
            updated_html = replace(html_content, r"const EMBEDDED_DATA = `[\s\S]*?`;" => new_data_js)
            write(html_path, updated_html)
            println("Updated $html_path with latest data.")
        end
    catch e
        println("Note: Could not auto-update HTML: $e")
    end
end

run_full_benchmark()
