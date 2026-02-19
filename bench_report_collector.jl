#!/usr/bin/env julia
"""
bench_report_collector.jl — High-Precision Performance Data Gatherer & Reporter.
Syncs results to CSV and auto-integrates into benchmark_results.html.
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

function run_report_benchmark()
    # Configuration: (n, m, trials)
    test_cases = [
        (1000, 5000, 100),
        (5000, 25000, 100),
        (10000, 50000, 50),
        (25000, 125000, 30),
        (50000, 250000, 20),
        (75000, 375000, 15),
        (100000, 500000, 15),
        (150000, 750000, 10),
        (200000, 1000000, 8),
        (250000, 1250000, 10),
        (350000, 1750000, 6),
        (500000, 2500000, 5),
        (750000, 3750000, 3),
        (1000000, 5000000, 2)
    ]

    println("DMMSY Synchronized Performance Reporter (v5.6)")
    
    # Correctness checks
    if !DMMSYSSSP.verify_correctness()
        println("Optimized Correctness check FAILED. Aborting.")
        return
    end
    
    println("===========================================================================")
    @printf("%-10s %-10s %-12s %-12s %-12s %-8s\n", "n", "m", "Dijkstra", "DMMSY Opt", "DMMSY Res", "Spd(Opt)")
    println("---------------------------------------------------------------------------")

    output_file = "benchmark_data.csv"
    
    open(output_file, "w") do f
        println(f, "n,m,dijkstra,dmmsy_opt,dmmsy_res")
        
        for (n, m, trials) in test_cases
            g = random_graph(n, m, 100.0)
            source = 1

            # Warm-up
            for _ in 1:2
                DijkstraModule.dijkstra_ref(g, source)
                DMMSYSSSP.ssp_duan(g, source)
                if n <= 50000; DMMSYResearch.ssp_duan_research(g, source) end
            end
            
            # Research version restricted to manageable sizes
            t_res_val = 0.0
            if n <= 100000
                t_res = @elapsed for _ in 1:trials DMMSYResearch.ssp_duan_research(g, source) end
                t_res_val = (t_res / trials) * 1000
            end

            # Execution for Optimized versions
            t_dij = @elapsed for _ in 1:trials DijkstraModule.dijkstra_ref(g, source) end
            t_opt = @elapsed for _ in 1:trials DMMSYSSSP.ssp_duan(g, source) end

            avg_dij = (t_dij / trials) * 1000
            avg_opt = (t_opt / trials) * 1000
            spd = avg_dij / avg_opt

            res_str = n <= 100000 ? @sprintf("%-12.4f", t_res_val) : "Skipped"
            @printf("%-10d %-10d %-12.4f %-12.4f %-12s %-8.2fx\n", n, m, avg_dij, avg_opt, res_str, spd)
            @printf(f, "%d,%d,%.6f,%.6f,%.6f\n", n, m, avg_dij, avg_opt, t_res_val)
        end
    end

    println("===========================================================================")
    println("Full benchmark complete. Results saved to $output_file")
    
    # Auto-embed into HTML for visualization
    try
        csv_data = read(output_file, String)
        html_path = "benchmark_results.html"
        if isfile(html_path)
            html_content = read(html_path, String)
            new_data_js = "const EMBEDDED_DATA = `$(csv_data)`;"
            updated_html = replace(html_content, r"const EMBEDDED_DATA = `[\s\S]*?`;" => new_data_js)
            write(html_path, updated_html)
            println("Updated $html_path with latest data.")
        end
    catch e
        println("Note: Could not auto-update HTML: $e")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_report_benchmark()
end
