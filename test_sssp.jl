if !@isdefined(Common) include("Common.jl") end
if !@isdefined(CSRGraphModule) include("CSRGraph.jl") end
if !@isdefined(DijkstraModule) include("Dijkstra.jl") end
if !@isdefined(DMMSYSSSP) include("DMMSY-SSSP.jl") end
if !@isdefined(DMMSYResearch) include("DMMSY_Research.jl") end

using .Common
using .CSRGraphModule
using .DijkstraModule
using .DMMSYSSSP
using .DMMSYResearch
using Test
using Random
Random.seed!(1234)

function run_comprehensive_tests()
    @testset "SSSP Correctness Tests" begin
        
        @testset "Topology: Simple Diamond" begin
            g = CSRGraph(4, [(1,2,2.0), (1,3,3.0), (2,4,1.0), (3,4,1.0)])
            d_ref, _ = dijkstra_ref(g, 1)
            
            d_opt, _ = ssp_duan(g, 1)
            @test d_opt ≈ d_ref
            
            d_res, _ = ssp_duan_research(g, 1)
            @test d_res ≈ d_ref
        end

        @testset "Topology: Cycle with shortcut" begin
            g = CSRGraph(3, [(1,2,1.0), (2,3,1.0), (3,2,0.1)])
            d_ref, _ = dijkstra_ref(g, 1)
            
            d_opt, _ = ssp_duan(g, 1)
            @test d_opt ≈ d_ref
            
            d_res, _ = ssp_duan_research(g, 1)
            @test d_res ≈ d_ref
        end

        @testset "Topology: Large Random Graph" begin
            # Using stable seed for reproducibility if needed, but rand is fine for general check
            n, m = 1000, 5000
            g = random_graph(n, m, 100.0)
            source = 1
            d_ref, _ = dijkstra_ref(g, source)
            
            d_opt, _ = ssp_duan(g, source)
            @test d_opt ≈ d_ref
            
            d_res, _ = ssp_duan_research(g, source)
            @test d_res ≈ d_ref
        end

        @testset "Topology: Disconnected Components" begin
            g = CSRGraph(5, [(1,2,1.0), (2,3,1.0), (4,5,1.0)])
            d_ref, _ = dijkstra_ref(g, 1)
            
            d_opt, _ = ssp_duan(g, 1)
            @test d_opt ≈ d_ref
            @test d_opt[4] == typemax(Float64)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_comprehensive_tests()
end
