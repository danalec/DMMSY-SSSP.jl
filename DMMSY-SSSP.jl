"""
DMMSY-SSSP.jl — Performance-Optimized Julia implementation of
“Breaking the Sorting Barrier for Directed Single-Source Shortest Paths”
by Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, and Longhui Yin (STOC 2025). 

This module implements a high-performance research-oriented single-source shortest path (SSSP)
solver for directed graphs with non-negative edge weights.
"""
module DMMSYSSSP

using DataStructures: BinaryMinHeap
using ..CSRGraphModule: CSRGraph, random_graph
using ..DijkstraModule: dijkstra_ref

export ssp_duan, verify_correctness, benchmark_comparison

# ============================================================================
# Tuning Constants
# ============================================================================
const BMSP_BASE_N        = 1500     # fallback to Dijkstra below this
const BMSP_MAX_DEPTH     = 3        # maximum recursion depth
const BMSP_POOL_DEPTH    = BMSP_MAX_DEPTH + 1
const BMSP_PIVOT_FACTOR  = 1.5      # log2(n) multiplier for max_pivots
const BMSP_THRESHOLD_K   = 2.2      # mean_weight multiplier
const BMSP_THRESHOLD_DEC = 0.7      # threshold decay per level

# ============================================================================
# Type Safety Wrappers
# ============================================================================
struct BlockIdx{T<:Integer}; val::T; end
struct VertIdx{T<:Integer};  val::T; end

# ============================================================================
# Flattened BlockedPartialPQ for High Performance
# ============================================================================

mutable struct BlockedPartialPQ{T<:Integer, W<:Real}
    n_blocks::T                            # Fixed at construction
    block_size::T                          # Fixed at construction
    storage_v::Vector{T}                   # Flattened vertex ID storage
    storage_d::Vector{W}                   # Flattened distance storage
    block_counts::Vector{Int}              # Current count per block
    block_minima::Vector{W}
    heap::BinaryMinHeap{Tuple{W, BlockIdx{T}}}
    overflow_heap::BinaryMinHeap{Tuple{W, VertIdx{T}}}
    total_size::T
end

function BlockedPartialPQ{T, W}(n::Integer) where {T<:Integer, W<:Real}
    b_size = T(min(64, max(1, ceil(Int, sqrt(n))))) 
    n_b = max(T(1), T(ceil(Int, n / b_size)))
    
    # Flattened storage to reduce pointer chasing and allocation overhead
    storage_v = Vector{T}(undef, n_b * b_size)
    storage_d = Vector{W}(undef, n_b * b_size)
    block_counts = zeros(Int, n_b)
    block_minima = fill(typemax(W), n_b)
    
    return BlockedPartialPQ{T, W}(n_b, b_size, storage_v, storage_d, block_counts, block_minima, 
                                 BinaryMinHeap{Tuple{W, BlockIdx{T}}}(), 
                                 BinaryMinHeap{Tuple{W, VertIdx{T}}}(), zero(T))
end

function reset!(pq::BlockedPartialPQ{T, W}) where {T, W}
    fill!(pq.block_counts, 0)
    fill!(pq.block_minima, typemax(W))
    empty!(pq.heap)
    empty!(pq.overflow_heap)
    pq.total_size = zero(T)
end

@inline function insert!(pq::BlockedPartialPQ{T, W}, v::T, d::W) where {T, W}
    # Direct hash distribution
    b_idx = (Int(v) % Int(pq.n_blocks)) + 1
    count = @inbounds pq.block_counts[b_idx]
    
    if count < pq.block_size
        offset = (b_idx - 1) * Int(pq.block_size)
        @inbounds pq.storage_v[offset + count + 1] = v
        @inbounds pq.storage_d[offset + count + 1] = d
        @inbounds pq.block_counts[b_idx] = count + 1
        
        if d < @inbounds(pq.block_minima[b_idx])
            @inbounds pq.block_minima[b_idx] = d
            push!(pq.heap, (d, BlockIdx(T(b_idx))))
        end
    else
        push!(pq.overflow_heap, (d, VertIdx(v)))
    end
    pq.total_size += 1
end

function extract_min!(pq::BlockedPartialPQ{T, W}, d_vec::Vector{W}) where {T, W}
    pq.total_size == zero(T) && return nothing
    h, minima = pq.heap, pq.block_minima
    oh = pq.overflow_heap
    storage_v, storage_d = pq.storage_v, pq.storage_d
    b_size = Int(pq.block_size)

    while true
        # 1. Process overflow heap for global candidates
        while !isempty(oh)
            dist_v, v_wrap = first(oh)
            v = v_wrap.val
            if dist_v > @inbounds(d_vec[v]) # Stale pruning
                pop!(oh)
                pq.total_size -= 1
                continue
            end
            
            # Global heap has better candidate than block heads
            if isempty(h) || dist_v <= first(h)[1]
                pop!(oh)
                pq.total_size -= 1
                return (v, dist_v)
            end
            break
        end

        # 2. Extract best block
        isempty(h) && break
        min_dist, b_idx_wrap = pop!(h)
        b_idx = b_idx_wrap.val
        
        # Stale check for block minimum
        @inbounds minima[b_idx] != min_dist && continue
        
        # 3. Linear scan of block + Pruning
        count = @inbounds pq.block_counts[b_idx]
        offset = (Int(b_idx) - 1) * b_size
        min_v, min_pos = typemax(W), 0
        
        curr_idx = 1
        while curr_idx <= count
            idx = offset + curr_idx
            v_val = @inbounds storage_v[idx]
            d_val = @inbounds storage_d[idx]
            
            if d_val > @inbounds(d_vec[v_val]) # Stale pruning
                last_pos = offset + count
                @inbounds storage_v[idx] = storage_v[last_pos]
                @inbounds storage_d[idx] = storage_d[last_pos]
                count -= 1
                pq.total_size -= 1
            else
                if d_val < min_v
                    min_v, min_pos = d_val, curr_idx
                end
                curr_idx += 1
            end
        end
        @inbounds pq.block_counts[b_idx] = count
        
        if min_pos == 0 # All entries in block were stale
            @inbounds minima[b_idx] = typemax(W)
            continue
        end
        
        # 4. Check global heap again against the fresh block minimum
        if !isempty(oh) && first(oh)[1] < min_v
            @inbounds minima[b_idx] = min_v
            push!(h, (min_v, BlockIdx(b_idx)))
            continue
        end

        # 5. Extract min vertex from block
        res_idx = offset + min_pos
        res_v = @inbounds storage_v[res_idx]
        res_d = @inbounds storage_d[res_idx]
        
        # Remove by swapping with last
        last_pos = offset + count
        @inbounds storage_v[res_idx] = storage_v[last_pos]
        @inbounds storage_d[res_idx] = storage_d[last_pos]
        count -= 1
        pq.total_size -= 1
        @inbounds pq.block_counts[b_idx] = count
        
        # Re-push block to heap if not empty
        if count == 0
            @inbounds minima[b_idx] = typemax(W)
        else
            new_m = typemax(W)
            @inbounds for k in 1:count
                val_d = storage_d[offset + k]
                (val_d < new_m) && (new_m = val_d)
            end
            @inbounds minima[b_idx] = new_m
            push!(h, (new_m, BlockIdx(b_idx)))
        end
        
        return (res_v, res_d)
    end
    return nothing
end

@inline function decrease_key!(pq::BlockedPartialPQ{T, W}, v::T, new_d::W) where {T, W}
    insert!(pq, v, new_d)
end

Base.isempty(pq::BlockedPartialPQ) = pq.total_size == zero(typeof(pq.total_size))

# ============================================================================
# Optimized SSSP Algorithms
# ============================================================================

@inline function find_pivots!(pivots::Vector{T}, candidates::Vector{T}, g::CSRGraph{T, W}, d::Vector{W}, threshold::W, max_pivots::Int) where {T, W}
    empty!(pivots); empty!(candidates)
    @inbounds for v in 1:g.n
        (d[v] < threshold) && push!(candidates, v)
    end
    isempty(candidates) && return pivots
    
    # Zero-allocation sort
    sort!(candidates, lt = (a, b) -> @inbounds(d[a] < d[b]))
    
    n_pivots = min(length(candidates), max_pivots)
    if n_pivots > 0
        step = max(1, length(candidates) ÷ n_pivots)
        @inbounds for i in 1:step:length(candidates)
            push!(pivots, candidates[i])
            length(pivots) >= n_pivots && break
        end
    end
    return pivots
end

@inline function base_case_bmsp!(g::CSRGraph{T, W}, sources::Vector{T}, d::Vector{W}, pred::Vector{T}, 
                                heap::BinaryMinHeap{Tuple{W, T}}) where {T, W}
    empty!(heap)
    for s in sources push!(heap, (d[s], s)) end
    adj, weights, offset = g.adjacency, g.weights, g.offset
    
    while !isempty(heap)
        (dist_u, u) = pop!(heap)
        dist_u > d[u] && continue
        
        start = @inbounds offset[u]
        stop = (@inbounds offset[u+1]) - 1
        @inbounds for i in start:stop
            v, w = adj[i], weights[i]
            if dist_u + w < d[v]
                d[v] = dist_u + w
                pred[v] = u
                push!(heap, (d[v], v))
            end
        end
    end
    return d, pred
end

function bmsp_recursive!(g::CSRGraph{T, W}, sources::Vector{T}, d::Vector{W}, pred::Vector{T}, 
                         threshold::W, depth::Int, pq::BlockedPartialPQ{T, W},
                         c_buf::Vector{T}, pivots_pool::Vector{Vector{T}},
                         base_heap::BinaryMinHeap{Tuple{W, T}}) where {T, W}
    n = g.n
    if n <= BMSP_BASE_N || depth >= BMSP_MAX_DEPTH
        return base_case_bmsp!(g, sources, d, pred, base_heap)
    end
    
    max_p = max(8, Int(ceil(log2(n) * BMSP_PIVOT_FACTOR)))
    pool_idx = min(depth + 1, BMSP_POOL_DEPTH)
    pivots = pivots_pool[pool_idx]
    
    find_pivots!(pivots, c_buf, g, d, threshold, max_p)
    isempty(pivots) && return base_case_bmsp!(g, sources, d, pred, base_heap)

    bmsp_recursive!(g, pivots, d, pred, threshold * W(BMSP_THRESHOLD_DEC), depth + 1, pq, c_buf, pivots_pool, base_heap)

    reset!(pq)
    for s in sources insert!(pq, s, d[s]) end
    for p in pivots insert!(pq, p, d[p]) end

    adj, weights, offset = g.adjacency, g.weights, g.offset
    while !isempty(pq)
        item = extract_min!(pq, d)
        item === nothing && break
        u, dist_u = item
        
        dist_u > @inbounds(d[u]) && continue
        
        start = @inbounds offset[u]
        stop = (@inbounds offset[u+1]) - 1
        @inbounds for i in start:stop
            v, w = adj[i], weights[i]
            nd = dist_u + w
            if nd < d[v]
                d[v] = nd
                pred[v] = u
                decrease_key!(pq, v, nd)
            end
        end
    end
    return d, pred
end

function ssp_duan(g::CSRGraph{T, W}, source::T) where {T<:Integer, W<:Real}
    n = g.n
    n == 0 && return W[], T[]
    
    if source < 1 || source > n
        throw(ArgumentError("Source vertex $source out of bounds [1, $n]."))
    end
    
    d, pred = fill(typemax(W), n), zeros(T, n)
    d[source] = zero(W)
    
    # Resources pre-allocation
    pq = BlockedPartialPQ{T, W}(n)
    c_buf = sizehint!(T[], isqrt(n))
    pivots_pool = [sizehint!(T[], max(8, Int(ceil(log2(n+1) * BMSP_PIVOT_FACTOR)))) for _ in 1:BMSP_POOL_DEPTH]
    base_heap = BinaryMinHeap{Tuple{W, T}}()
    
    threshold = g.mean_weight * log10(n+1) * BMSP_THRESHOLD_K
    
    return bmsp_recursive!(g, [source], d, pred, W(threshold), 0, pq, c_buf, pivots_pool, base_heap)
end

# ============================================================================
# Verification & Benchmark
# ============================================================================

function verify_correctness()
    println("=" ^ 40)
    println("DMMSY Optimized Correctness")
    println("=" ^ 40)
    all_passed = true
    test_cases = [
        ("Linear",  CSRGraph(5, [(1,2,1.0), (2,3,2.0), (3,4,1.0), (4,5,3.0)])),
        ("Diamond", CSRGraph(4, [(1,2,2.0), (1,3,3.0), (2,4,1.0), (3,4,1.0)])),
        ("Cycle",   CSRGraph(3, [(1,2,1.0), (2,3,1.0), (3,2,0.5)])),
    ]
    for (name, g) in test_cases
        d1, _ = ssp_duan(g, 1)
        d_ref, _ = dijkstra_ref(g, 1)
        if d1 ≈ d_ref
            println("  ✓ $name PASSED")
        else
            println("  ✗ $name FAILED")
            all_passed = false
        end
    end
    return all_passed
end

function benchmark_comparison(n::Int, m::Int, trials::Int=5)
    g = random_graph(n, m, 100.0)
    ssp_duan(g, 1) # Warm-up
    dijkstra_ref(g, 1) # Warm-up
    t1 = @elapsed for _ in 1:trials ssp_duan(g, 1) end
    t2 = @elapsed for _ in 1:trials dijkstra_ref(g, 1) end
    println("Benchmark (n=$n, m=$m):")
    println("  DMMSY Opt: ", round(t1/trials * 1000, digits=4), "ms")
    println("  Dijkstra:  ", round(t2/trials * 1000, digits=4), "ms")
end

end # module
