
#!/usr/bin/env julia
"""
bench_report_collector.jl — High-Precision Performance Data Gatherer & Reporter.
Syncs results to CSV and auto-integrates into benchmark_results.html.
"""

using Random
Random.seed!(1234)

include("Common.jl")
include("CSRGraph.jl")
include("Dijkstra.jl")
include("DMMSY-SSSP.jl")
include("DMMSY_Research.jl")

using .Common
using .CSRGraphModule
using .DijkstraModule
using .DMMSYSSSP
using .DMMSYResearch
using Printf
include("test_sssp.jl")

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
        (500000, 2500000, 20),
        (750000, 3750000, 15),
        (1000000, 5000000, 10)
    ]

    println("DMMSY Synchronized Performance Reporter (v5.6)")

    # Correctness checks via external test suite
    try
        run_comprehensive_tests()
        println("All correctness checks PASSED.")
    catch e
        println("Correctness check FAILED: $e. Aborting.")
        return
    end

    println("===========================================================================")
    @printf("%-10s %-10s %-12s %-12s %-12s %-8s\n", "n", "m", "Dijkstra", "DMMSY Res", "DMMSY Opt", "Spd(Opt)")
    println("---------------------------------------------------------------------------")

    output_file = "benchmark_data.csv"

    open(output_file, "w") do f
        println(f, "n,m,dijkstra,dmmsy_res,dmmsy_opt")

        for (n, m, trials) in test_cases
            g = random_graph(n, m, 100.0)
            source = 1

            # Warm-up for all implementations (fair comparison)
            for _ in 1:min(5, trials ÷ 2)
                DijkstraModule.dijkstra_ref(g, source)
                DMMSYSSSP.ssp_duan(g, source)
                DMMSYResearch.ssp_duan_research(g, source)
            end

            # Execution for all implementations with consistent trials
            GC.gc()
            GC.enable(false)
            t_dij = @elapsed for _ in 1:trials DijkstraModule.dijkstra_ref(g, source) end
            GC.enable(true)



            GC.gc()
            GC.enable(false)
            t_opt = @elapsed for _ in 1:trials DMMSYSSSP.ssp_duan(g, source) end
            GC.enable(true)

            GC.gc()
            GC.enable(false)
            t_res = @elapsed for _ in 1:trials DMMSYResearch.ssp_duan_research(g, source) end
            GC.enable(true)

            avg_dij = (t_dij / trials) * 1000
            avg_opt = (t_opt / trials) * 1000
            avg_res = (t_res / trials) * 1000
            spd = avg_dij / avg_opt

            @printf("%-10d %-10d %-12.4f %-12.4f %-12.4f %-8.2fx\n", n, m, avg_dij, avg_res, avg_opt, spd)
            @printf(f, "%d,%d,%.6f,%.6f,%.6f\n", n, m, avg_dij, avg_res, avg_opt)
        end
    end

    println("===========================================================================")
    println("Full benchmark complete. Results saved to $output_file")

    # Write JS file so HTML can load the data without a local server (CORS workaround)
    try
        csv_data = read(output_file, String)
        open("benchmark_data.js", "w") do f
            write(f, "const BENCHMARK_CSV_DATA = \`\n$(csv_data)\`;\n")
        end
        println("Exported JavaScript data binding to benchmark_data.js")
    catch e
        println("Note: Could not create JS data file: $e")
    end

    println("Data is written. You can now open `benchmark_results.html` directly in your browser without a local server.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_report_benchmark()
end