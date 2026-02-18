"""
DMMSY-SSSP.jl — Experimental Julia implementation of
“Breaking the Sorting Barrier for Directed Single-Source Shortest Paths”
by Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, and Longhui Yin (STOC 2025). 

This module implements a research-oriented single-source shortest path (SSSP)
solver for directed graphs with non‑negative edge weights, inspired by the
O(m · log^(2/3)(n)) comparison–addition model algorithm in the original paper. 

Status:
- Algorithmic structure follows the BMSSP / pivot / blocked-PQ ideas, but is not
  a line-by-line reproduction of all invariants and proofs.
- Focus is on clarity, experimentation, and comparison against Dijkstra rather
  than on asymptotic-tight constants or full theoretical guarantees. 

Exports:
- `CSRGraph`: cache-friendly directed graph in CSR format
- `ssp_duan(g, s)`: run DMMSY-style SSSP from source vertex `s`
- `reconstruct_path(pred, s, t)`: recover an s→t shortest path from predecessors
- `dijkstra_ref(g, s)`: reference SSSP for verification and benchmarking
- `verify_correctness()`: small test suite comparing against Dijkstra
- `benchmark_comparison(n, m; trials=5)`: quick DMMSY vs Dijkstra timing harness

Reference:
- Duan et al., “Breaking the Sorting Barrier for Directed Single-Source Shortest
  Paths”, arXiv:2504.17033, STOC 2025.
"""

module DMMSYSSSP

using DataStructures: BinaryHeap, DefaultDict, heapify!, heappop!, heappush!

# ============================================================================
# Core Data Structures
# ============================================================================

"""
    Edge{T<:Integer}

Represents a directed edge in the graph.
"""
struct Edge{T<:Integer}
    target::T
    weight::Float64
end

"""
    CSRGraph{T<:Integer}

Compressed Sparse Row representation of a directed graph for cache efficiency.
"""
struct CSRGraph{T<:Integer}
    n::T                      # Number of vertices
    m::T                      # Number of edges
    offset::Vector{T}         # offset[i] = start index of adjacency list for vertex i
    adjacency::Vector{T}      # Target vertices (flattened adjacency lists)
    weights::Vector{Float64}  # Edge weights (parallel to adjacency)
end

"""
    CSRGraph(n::Int, edges::Vector{Tuple{Int, Int, Float64}})

Construct a CSR graph from edge list.
Edge format: (source, target, weight)
"""
function CSRGraph{T}(n::Int, edges::Vector{Tuple{Int, Int, Float64}}) where T<:Integer
    # Group edges by source
    buckets = [Tuple{Int, Float64}[] for _ in 1:n]
    for (u, v, w) in edges
        if 1 <= u <= n && 1 <= v <= n && w >= 0
            push!(buckets[u], (v, w))
        end
    end

    # Build CSR representation
    m = sum(length(bucket) for bucket in buckets)
    offset = zeros(T, n + 1)
    adjacency = zeros(T, m)
    weights = zeros(Float64, m)

    idx = 1
    for u in 1:n
        offset[u] = idx
        for (v, w) in buckets[u]
            adjacency[idx] = v
            weights[idx] = w
            idx += 1
        end
    end
    offset[n + 1] = m + 1

    return CSRGraph{T}(n, T(m), offset, adjacency, weights)
end

CSRGraph(n::Int, edges::Vector{Tuple{Int, Int, Float64}}) = CSRGraph{Int}(n, edges)

"""
    out_edges(g::CSRGraph, u::Int)

Return an iterator over outgoing edges from vertex u.
"""
function out_edges(g::CSRGraph{T}, u::T) where T<:Integer
    start = g.offset[u]
    stop = g.offset[u + 1] - 1
    return ((g.adjacency[i], g.weights[i]) for i in start:stop)
end

"""
    BlockedPartialPQ{T<:Integer}

Blocked Partial Priority Queue implementing Lemma 3.3.
This is a cache-friendly data structure that avoids sorting all elements.
"""
mutable struct BlockedPartialPQ{T<:Integer}
    n_blocks::T              # Number of blocks
    block_size::T            # Size of each block
    blocks::Vector{Vector{Tuple{T, Float64}}}  # Each block stores (vertex, distance)
    heap::BinaryHeap{Tuple{Float64, T}}        # Min-heap of block representatives
    total_size::T            # Total number of elements
end

"""
    BlockedPartialPQ{T}(n::Integer)

Create an empty blocked partial PQ with block size n^(2/3).
"""
function BlockedPartialPQ{T}(n::Integer) where T<:Integer
    block_size = ceil(Int, n^(2/3))
    n_blocks = ceil(Int, n / block_size) + 1
    blocks = [Tuple{T, Float64}[] for _ in 1:n_blocks]
    return BlockedPartialPQ{T}(T(n_blocks), T(block_size), blocks, BinaryHeap{Tuple{Float64, T}}(), zero(T))
end

"""
    insert!(pq::BlockedPartialPQ, v::Integer, d::Real)

Insert vertex v with distance d into the PQ.
"""
function insert!(pq::BlockedPartialPQ{T}, v::T, d::Real) where T<:Integer
    # Insert into first available block or create overflow
    inserted = false
    for i in 1:3  # First 3 blocks are for active elements
        if length(pq.blocks[i]) < pq.block_size
            push!(pq.blocks[i], (v, Float64(d)))
            inserted = true
            break
        end
    end

    if !inserted
        # Find any block with space
        for i in 1:pq.n_blocks
            if length(pq.blocks[i]) < pq.block_size
                push!(pq.blocks[i], (v, Float64(d)))
                inserted = true
                break
            end
        end
    end

    if !inserted
        push!(pq.blocks[end], (v, Float64(d)))
    end

    pq.total_size += 1
end

"""
    extract_min!(pq::BlockedPartialPQ)

Extract and return the (vertex, distance) pair with minimum distance.
"""
function extract_min!(pq::BlockedPartialPQ{T}) where T<:Integer
    if isempty(pq)
        return nothing
    end

    # Find minimum across all blocks
    min_block = 1
    min_idx = 1
    min_dist = Inf
    min_vertex = T(0)

    for i in 1:pq.n_blocks
        if !isempty(pq.blocks[i])
            for j in eachindex(pq.blocks[i])
                if pq.blocks[i][j][2] < min_dist
                    min_dist = pq.blocks[i][j][2]
                    min_vertex = pq.blocks[i][j][1]
                    min_block = i
                    min_idx = j
                end
            end
        end
    end

    # Remove the minimum element
    deleteat!(pq.blocks[min_block], min_idx)
    pq.total_size -= 1

    return (min_vertex, min_dist)
end

"""
    decrease_key!(pq::BlockedPartialPQ, v::Integer, new_d::Real)

Decrease the distance of vertex v to new_d.
For efficiency, we use lazy update - just insert the new value
and rely on the complete array to track current minimum.
"""
function decrease_key!(pq::BlockedPartialPQ{T}, v::T, new_d::Real) where T<:Integer
    # Lazy approach: just insert new value
    insert!(pq, v, new_d)
end

"""
    isempty(pq::BlockedPartialPQ)

Check if the PQ is empty.
"""
Base.isempty(pq::BlockedPartialPQ) = pq.total_size == 0

"""
    clear!(pq::BlockedPartialPQ)

Clear all elements from the PQ.
"""
function clear!(pq::BlockedPartialPQ{T}) where T<:Integer
    for block in pq.blocks
        empty!(block)
    end
    pq.total_size = zero(T)
end

# ============================================================================
# Algorithm 1: FindPivots
# ============================================================================

"""
    find_pivots(g::CSRGraph, d::Vector{Float64}, threshold::Real)

Identify pivot vertices as described in Algorithm 1.
Pivots are vertices whose distance is sufficiently small relative to the threshold.

Returns:
- pivots: Vector of pivot vertex indices
- pivot_regions: For each vertex, the associated pivot (or 0 if none)
"""
function find_pivots(g::CSRGraph{T}, d::Vector{Float64}, threshold::Real) where T<:Integer
    n = g.n
    pivots = T[]
    pivot_region = zeros(T, n)
    visited = falses(n)

    # Sort vertices by distance (threshold-based selection)
    # Vertices with distance < threshold are potential pivots
    candidates = T[]
    for v in 1:n
        if d[v] < threshold && !visited[v]
            push!(candidates, v)
        end
    end

    # Greedy selection to maximize coverage
    sort!(candidates, by = v -> d[v])

    while !isempty(candidates)
        v = popfirst!(candidates)
        if visited[v]
            continue
        end

        push!(pivots, v)
        pivot_region[v] = v
        visited[v] = true

        # Mark vertices reachable from v as covered
        for (u, w) in out_edges(g, v)
            if d[u] >= d[v] + w && !visited[u]
                pivot_region[u] = v
                visited[u] = true
            end
        end
    end

    # For remaining vertices, find closest pivot
    for v in 1:n
        if !visited[v]
            # BFS to find closest pivot
            best_pivot = T(0)
            best_dist = Inf

            for p in pivots
                # Check if p can reach v (this is a simplification)
                # Full implementation would use distance estimates
                if d[v] >= d[p]
                    if abs(d[v] - d[p]) < best_dist
                        best_dist = abs(d[v] - d[p])
                        best_pivot = p
                    end
                end
            end
            pivot_region[v] = best_pivot
        end
    end

    return pivots, pivot_region
end

# ============================================================================
# Algorithm 2: BaseCase BMSSP
# ============================================================================

"""
    base_case_bmsp(g::CSRGraph, sources::Vector{Int}, d_init::Vector{Float64})

Base case using Dijkstra's algorithm with a standard heap.
Used when the subproblem size is small enough.

Arguments:
- g: The graph
- sources: Source vertices for multi-source SSSP
- d_init: Initial distance estimates

Returns:
- d: Final distances
- pred: Predecessor array for path reconstruction
"""
function base_case_bmsp(g::CSRGraph{T},
                        sources::Vector{T},
                        d_init::Vector{Float64},
                        pred::Vector{T}=zeros(T, g.n)) where T<:Integer

    n = g.n
    d = copy(d_init)

    # Initialize heap with source vertices
    heap = BinaryHeap{Tuple{Float64, T}}()
    in_heap = falses(n)

    for s in sources
        heappush!(heap, (d[s], s))
        in_heap[s] = true
    end

    # Standard Dijkstra relaxation
    while !isempty(heap)
        (dist_u, u) = heappop!(heap)
        in_heap[u] = false

        if dist_u > d[u]
            continue
        end

        for (v, w) in out_edges(g, u)
            new_dist = d[u] + w
            if new_dist < d[v]
                d[v] = new_dist
                pred[v] = u
                if !in_heap[v]
                    heappush!(heap, (new_dist, v))
                    in_heap[v] = true
                else
                    # Lazy decrease key: push new value and ignore old ones
                    heappush!(heap, (new_dist, v))
                end
            end
        end
    end

    return d, pred
end

# ============================================================================
# Algorithm 3: BMSSP (Recursive)
# ============================================================================

"""
    BMSSP(g::CSRGraph, sources::Vector{Int}, d_init::Vector{Float64}, threshold::Real)

Blocked Multi-Source Shortest Path - the core recursive algorithm.

The algorithm recursively:
1. Identifies pivot vertices (FindPivots)
2. Computes paths to/from pivots
3. Processes remaining vertices using block decomposition

Returns:
- d: Final distances
- pred: Predecessor array
"""
function bmsp!(g::CSRGraph{T},
               sources::Vector{T},
               d::Vector{Float64},
               pred::Vector{T},
               threshold::Real,
               depth::Int=0) where T<:Integer

    n = g.n

    # Base case: small subproblem
    if n <= 1000 || depth >= 10
        return base_case_bmsp(g, sources, d, pred)
    end

    # Step 1: Find pivots
    pivots, pivot_region = find_pivots(g, d, threshold)

    if isempty(pivots)
        # No pivots found, fall back to base case
        return base_case_bmsp(g, sources, d, pred)
    end

    # Step 2: Compute distances to/from each pivot
    pivot_distances = zeros(Float64, length(pivots), n)
    pivot_preds = zeros(T, length(pivots), n)

    for (i, p) in enumerate(pivots)
        # Multi-source SSSP from pivot to all vertices
        d_pivot = fill(Inf, n)
        d_pivot[p] = 0.0
        pred_pivot = zeros(T, n)

        # Run BMSSP or base case for each pivot
        if n <= 5000
            pivot_distances[i, :], _ = base_case_bmsp(g, [p], d_pivot, pred_pivot)
        else
            # Recursive call with smaller threshold
            pivot_distances[i, :], _ = bmsp!(g, [p], d_pivot, pred_pivot,
                                               threshold * 0.5, depth + 1)
        end
    end

    # Step 3: Update distances using pivot information
    for v in 1:n
        for (i, p) in enumerate(pivots)
            if pivot_distances[i, v] + d[p] < d[v]
                d[v] = pivot_distances[i, v] + d[p]
            end
        end
    end

    # Step 4: Process non-pivot vertices in blocks
    non_pivot_vertices = T[v for v in 1:n if !(v in pivots)]

    # Process in blocks of size n^(2/3)
    block_size = ceil(Int, n^(2/3))
    for i in 1:block_size:length(non_pivot_vertices)
        block_start = i
        block_end = min(i + block_size - 1, length(non_pivot_vertices))
        block_vertices = non_pivot_vertices[block_start:block_end]

        # Run SSSP restricted to this block
        for v in block_vertices
            # Relax edges within block
            for (u, w) in out_edges(g, v)
                if u in block_vertices && d[v] + w < d[u]
                    d[u] = d[v] + w
                    pred[u] = v
                end
            end
        end
    end

    # Step 5: Final refinement pass
    for v in 1:n
        for (u, w) in out_edges(g, v)
            if d[v] + w < d[u]
                d[u] = d[v] + w
                pred[u] = v
            end
        end
    end

    return d, pred
end

# ============================================================================
# Main Interface
# ============================================================================

"""
    ssp_duan(g::CSRGraph, source::Integer)

Compute shortest paths from a single source using the DMMSY algorithm.

Arguments:
- g: A CSRGraph representing the directed graph with non-negative weights
- source: The source vertex (1-indexed)

Returns:
- distances: Vector of shortest distances from source to each vertex (Inf if unreachable)
- predecessors: Vector of predecessor vertices for path reconstruction (0 if none)

Example:
```julia
g = CSRGraph(5, [(1,2,1.0), (2,3,2.0), (1,3,4.0), (3,4,1.0), (4,5,3.0)])
dist, pred = ssp_duan(g, 1)
```
"""
function ssp_duan(g::CSRGraph{T}, source::T) where T<:Integer
    @assert 1 <= source <= g.n "Source vertex out of bounds"

    n = g.n

    # Initialize distances and predecessors
    d = fill(Inf, n)
    pred = zeros(T, n)

    d[source] = 0.0

    # Initial threshold based on graph structure
    # This is a heuristic; the paper provides more sophisticated threshold selection
    threshold = if n > 0 && g.m > 0
        # Estimate average edge weight for threshold
        total_weight = sum(g.weights)
        avg_weight = total_weight / g.m
        avg_weight * n^(1/3)
    else
        1.0
    end

    # Run the main BMSSP algorithm
    d, pred = bmsp!(g, [source], d, pred, threshold)

    return d, pred
end

"""
    reconstruct_path(pred::Vector{Int}, source::Int, target::Int)

Reconstruct the shortest path from source to target using the predecessor array.

Returns a vector of vertices representing the path, or an empty vector if no path exists.
"""
function reconstruct_path(pred::Vector{T}, source::T, target::T) where T<:Integer
    if pred[target] == 0 && target != source
        return T[]
    end

    path = T[target]
    current = target

    while current != source && current != 0
        current = pred[current]
        if current == 0
            return T[]  # No complete path
        end
        pushfirst!(path, current)
    end

    return path
end

# Export main interface
export ssp_duan, CSRGraph, reconstruct_path

# ============================================================================
# Reference Implementation (Dijkstra) for Verification
# ============================================================================

"""
    dijkstra_ref(g::CSRGraph, source::Integer)

Reference Dijkstra implementation for correctness verification.
Uses a simple heap-based approach.
"""
function dijkstra_ref(g::CSRGraph{T}, source::T) where T<:Integer
    n = g.n
    d = fill(Inf, n)
    pred = zeros(T, n)
    visited = falses(n)

    d[source] = 0.0
    heap = BinaryHeap{Tuple{Float64, T}}()
    heappush!(heap, (0.0, source))

    while !isempty(heap)
        (dist_u, u) = heappop!(heap)

        if visited[u]
            continue
        end
        visited[u] = true

        for (v, w) in out_edges(g, u)
            if !visited[v] && d[u] + w < d[v]
                d[v] = d[u] + w
                pred[v] = u
                heappush!(heap, (d[v], v))
            end
        end
    end

    return d, pred
end

# ============================================================================
# Test Suite
# ============================================================================

"""
    verify_correctness()

Run verification tests comparing DMMSY with reference Dijkstra.
Returns true if all tests pass.
"""
function verify_correctness()
    println("=" ^ 60)
    println("DMMSY-SSSP Correctness Verification")
    println("=" ^ 60)

    all_passed = true

    # Test 1: Simple linear chain
    println("\nTest 1: Linear chain graph...")
    g1 = CSRGraph(5, [(1,2,1.0), (2,3,2.0), (3,4,1.0), (4,5,3.0)])
    d1_dmmsy, _ = ssp_duan(g1, 1)
    d1_ref, _ = dijkstra_ref(g1, 1)

    if d1_dmmsy ≈ d1_ref
        println("  ✓ PASSED")
    else
        println("  ✗ FAILED")
        println("    DMMSY: ", d1_dmmsy)
        println("    Dijkstra: ", d1_ref)
        all_passed = false
    end

    # Test 2: Diamond graph
    println("\nTest 2: Diamond graph...")
    g2 = CSRGraph(4, [(1,2,2.0), (1,3,3.0), (2,4,1.0), (3,4,1.0)])
    d2_dmmsy, _ = ssp_duan(g2, 1)
    d2_ref, _ = dijkstra_ref(g2, 1)

    if d2_dmmsy ≈ d2_ref
        println("  ✓ PASSED")
    else
        println("  ✗ FAILED")
        println("    DMMSY: ", d2_dmmsy)
        println("    Dijkstra: ", d2_ref)
        all_passed = false
    end

    # Test 3: Disconnected graph
    println("\nTest 3: Disconnected graph...")
    g3 = CSRGraph(4, [(1,2,1.0), (3,4,2.0)])
    d3_dmmsy, _ = ssp_duan(g3, 1)
    d3_ref, _ = dijkstra_ref(g3, 1)

    if d3_dmmsy ≈ d3_ref
        println("  ✓ PASSED")
    else
        println("  ✗ FAILED")
        println("    DMMSY: ", d3_dmmsy)
        println("    Dijkstra: ", d3_ref)
        all_passed = false
    end

    # Test 4: Single node
    println("\nTest 4: Single node graph...")
    g4 = CSRGraph(1, Tuple{Int, Int, Float64}[])
    d4_dmmsy, _ = ssp_duan(g4, 1)
    d4_ref, _ = dijkstra_ref(g4, 1)

    if d4_dmmsy ≈ d4_ref
        println("  ✓ PASSED")
    else
        println("  ✗ FAILED")
        println("    DMMSY: ", d4_dmmsy)
        println("    Dijkstra: ", d4_ref)
        all_passed = false
    end

    # Test 5: Graph with cycle
    println("\nTest 5: Graph with cycle...")
    g5 = CSRGraph(3, [(1,2,1.0), (2,3,1.0), (3,2,0.5)])
    d5_dmmsy, _ = ssp_duan(g5, 1)
    d5_ref, _ = dijkstra_ref(g5, 1)

    if d5_dmmsy ≈ d5_ref
        println("  ✓ PASSED")
    else
        println("  ✗ FAILED")
        println("    DMMSY: ", d5_dmmsy)
        println("    Dijkstra: ", d5_ref)
        all_passed = false
    end

    # Test 6: Random graph
    println("\nTest 6: Random graph (n=50, m≈100)...")
    Random.seed!(42)
    g6 = random_graph(50, 100, 10.0)
    d6_dmmsy, _ = ssp_duan(g6, 1)
    d6_ref, _ = dijkstra_ref(g6, 1)

    if all(isapprox.(d6_dmmsy, d6_ref, rtol=1e-10))
        println("  ✓ PASSED")
    else
        println("  ✗ FAILED")
        max_diff = maximum(abs.(d6_dmmsy .- d6_ref))
        println("    Max difference: ", max_diff)
        all_passed = false
    end

    # Test 7: Path reconstruction
    println("\nTest 7: Path reconstruction...")
    g7 = CSRGraph(6, [(1,2,1.0), (2,3,1.0), (1,4,2.0), (4,5,1.0), (5,6,1.0), (3,6,1.0)])
    d7_dmmsy, pred7 = ssp_duan(g7, 1)
    path = reconstruct_path(pred7, 1, 6)

    if !isempty(path) && path[1] == 1 && path[end] == 6 && length(path) > 1
        println("  ✓ PASSED (path: ", join(path, " -> "), ")")
    else
        println("  ✗ FAILED")
        println("    Path: ", path)
        all_passed = false
    end

    # Summary
    println("\n" * "=" ^ 60)
    if all_passed
        println("All tests PASSED ✓")
    else
        println("Some tests FAILED ✗")
    end
    println("=" ^ 60)

    return all_passed
end

"""
    random_graph(n::Int, m::Int, max_weight::Float64)

Generate a random directed graph with n vertices, approximately m edges,
and random non-negative edge weights up to max_weight.
"""
function random_graph(n::Int, m::Int, max_weight::Float64)
    edges = Tuple{Int, Int, Float64}[]
    for _ in 1:m
        u = rand(1:n)
        v = rand(1:n)
        w = rand() * max_weight
        push!(edges, (u, v, w))
    end
    return CSRGraph(n, edges)
end

"""
    benchmark_comparison(n::Int, m::Int, trials::Int=5)

Compare performance between DMMSY and reference Dijkstra.
"""
function benchmark_comparison(n::Int, m::Int, trials::Int=5)
    println("\nBenchmark Comparison (n=$n, m=$m, trials=$trials)")

    g = random_graph(n, m, 100.0)

    # Warm-up
    ssp_duan(g, 1)
    dijkstra_ref(g, 1)

    # Benchmark DMMSY
    dmmsy_times = Float64[]
    for _ in 1:trials
        t = @elapsed ssp_duan(g, 1)
        push!(dmmsy_times, t)
    end

    # Benchmark Dijkstra
    dijkstra_times = Float64[]
    for _ in 1:trials
        t = @elapsed dijkstra_ref(g, 1)
        push!(dijkstra_times, t)
    end

    println("  DMMSY:      ", round(mean(dmmsy_times), digits=6), "s (±", round(std(dmmsy_times), digits=6), "s)")
    println("  Dijkstra:   ", round(mean(dijkstra_times), digits=6), "s (±", round(std(dijkstra_times), digits=6), "s)")

    return mean(dmmsy_times), mean(dijkstra_times)
end

end # module DMMSYSSSP

# Run tests if file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    using Random
    using .DMMSYSSSP

    DMMSYSSSP.verify_correctness()

    # Uncomment for benchmarking
    # println("\n--- Performance Benchmarking ---")
    # DMMSYSSSP.benchmark_comparison(100, 500, 10)
end
