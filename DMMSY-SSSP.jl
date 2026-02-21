"""
DMMSY-SSSP.jl — State-of-the-art SSSP.
Based on "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths"
Optimized implementation with bound filtering and efficient recursion.
"""
module DMMSYSSSP

using ..CSRGraphModule: CSRGraph
using ..DijkstraModule: dijkstra_ref
using ..Common: get_params, Edge, HeapNode, push_dec!, pop_min!, Fast4AryHeap

export ssp_duan, verify_correctness

# ----------------------------------------------------------------------------
# Optimized Workspace with pre-allocated buffers
# ----------------------------------------------------------------------------
mutable struct WorkSpace{T<:Integer, W<:Real}
    d::Vector{W}
    pr::Vector{T}
    h_nodes::Vector{HeapNode{T, W}}
    h_pos::Vector{T}
    dirty_h::Vector{T}
    dh_cnt::Int
    dirty_d::Vector{Int}
    ds_cnt::Int
    piv_bufs::Vector{Vector{T}}
end

function get_workspace(n::Int, ::Type{T}, ::Type{W}) where {T, W}
    tls = task_local_storage()
    key = (:dmmsy_v570_ws, T, W)
    ws = get!(tls, key) do
        WorkSpace{T, W}(
            fill(typemax(W), n),
            zeros(T, n),
            Vector{HeapNode{T, W}}(undef, n),
            zeros(T, n),
            zeros(T, n),
            0,
            zeros(Int, n),
            0,
            [Vector{T}(undef, max(4, Int(floor(log2(n)^(1/3))))) for _ in 1:max(2, Int(floor(log2(n)^(2/3)))) + 1]
        )
    end::WorkSpace{T, W}
    if length(ws.d) != n
        ws = WorkSpace{T, W}(
            fill(typemax(W), n),
            zeros(T, n),
            Vector{HeapNode{T, W}}(undef, n),
            zeros(T, n),
            zeros(T, n),
            0,
            zeros(Int, n),
            0,
            [Vector{T}(undef, max(4, Int(floor(log2(n)^(1/3))))) for _ in 1:max(2, Int(floor(log2(n)^(2/3)))) + 1]
        )
        tls[key] = ws
    end
    return ws
end

# ----------------------------------------------------------------------------
# Optimized BMSSP Implementation (Algorithm 3)
# ----------------------------------------------------------------------------
@inline function bmsp_rec!(g::CSRGraph{T,W}, d::Vector{W}, pr::Vector{T},
                          src_buf::Vector{T}, off_src::Int, len_src::Int,
                          B::W, dp::Int, ws::WorkSpace{T,W}) where {T,W}
    n = g.n
    k, t = get_params(n)
    off, edges = g.offset, g.edges

    # Base case: dp >= t or small src
    if dp >= t || len_src <= k
        # Selective reset of heap position map
        if ws.dh_cnt > 0
            @inbounds for i in 1:ws.dh_cnt
                ws.h_pos[Int(ws.dirty_h[i])] = zero(T)
            end
            ws.dh_cnt = 0
        end
        
        h = Fast4AryHeap{T, W}(ws.h_nodes, ws.h_pos, ws.dirty_h)
        sz, dcnt = 0, 0

        for i in 1:len_src
            @inbounds s = src_buf[off_src + i]
            sz, dcnt = push_dec!(h, sz, dcnt, s, d[Int(s)])
        end
        ws.dh_cnt = dcnt

        while sz > 0
            du, u, sz = pop_min!(h, sz)
            u_i = Int(u)
            @inbounds if du > d[u_i]; continue end

            @inbounds u_off = off[u_i]
            @inbounds u_end = off[u_i+1] - 1
            
            # Relaxation loop
            @fastmath for i in u_off:u_end
                @inbounds e = edges[i]
                @inbounds v, w = e.v, e.w
                v_i = Int(v)
                nd = du + w
                @inbounds if nd < d[v_i]
                    if d[v_i] == typemax(W)
                        ws.ds_cnt += 1
                        ws.dirty_d[ws.ds_cnt] = v_i
                    end
                    d[v_i], pr[v_i] = nd, u
                    sz, dcnt = push_dec!(h, sz, ws.dh_cnt, v, nd)
                    ws.dh_cnt = dcnt
                end
            end
        end
        return
    end

    # Select k pivots evenly from src
    np = min(len_src, k)
    pivots = ws.piv_bufs[dp + 2]
    step = max(1, len_src ÷ np)
    curr_np = 0
    bound = min(len_src, step * k)
    for i in 1:step:bound
        curr_np += 1
        @inbounds pivots[curr_np] = src_buf[off_src + i]
    end

    # Recursive call on pivots with tighter bound - zero allocation passing
    bmsp_rec!(g, d, pr, pivots, 0, curr_np, B * W(0.5), dp + 1, ws)

    # Main loop - process remaining src vertices with bound B
    # Selective reset of heap position map
    if ws.dh_cnt > 0
        @inbounds for i in 1:ws.dh_cnt
            ws.h_pos[Int(ws.dirty_h[i])] = zero(T)
        end
        ws.dh_cnt = 0
    end
    
    h = Fast4AryHeap{T, W}(ws.h_nodes, ws.h_pos, ws.dirty_h)
    sz, dcnt = 0, 0

    has_work = false
    @inbounds for i in 1:len_src
        s = src_buf[off_src + i]
        s_i = Int(s)
        dv = d[s_i]
        if dv < B
            sz, dcnt = push_dec!(h, sz, dcnt, s, dv)
            has_work = true
        end
    end
    ws.dh_cnt = dcnt

    if !has_work; return end

    while sz > 0
        du, u, sz = pop_min!(h, sz)
        u_i = Int(u)
        @inbounds if du > d[u_i]; continue end

        @inbounds u_off = off[u_i]
        @inbounds u_end = off[u_i+1] - 1
        @fastmath for i in u_off:u_end
            @inbounds e = edges[i]
            @inbounds v, w = e.v, e.w
            v_i = Int(v)
            nd = du + w
            @inbounds if nd < d[v_i]
                if d[v_i] == typemax(W)
                    ws.ds_cnt += 1
                    ws.dirty_d[ws.ds_cnt] = v_i
                end
                d[v_i], pr[v_i] = nd, u
                if nd < B
                    sz, dcnt = push_dec!(h, sz, ws.dh_cnt, v, nd)
                    ws.dh_cnt = dcnt
                end
            end
        end
    end
end

# ----------------------------------------------------------------------------
# Main Interface Function
# ----------------------------------------------------------------------------
function ssp_duan(g::CSRGraph{T, W}, src::T) where {T, W}
    n = g.n
    if n == 0
        return W[], T[]
    end

    ws = get_workspace(n, T, W)

    # Initialize distances - use selective reset
    if ws.ds_cnt > (n >> 2)
        fill!(ws.d, typemax(W))
        fill!(ws.pr, zero(T))
    else
        @inbounds for i in 1:ws.ds_cnt
            idx = ws.dirty_d[i]
            ws.d[idx] = typemax(W)
            ws.pr[idx] = zero(T)
        end
    end

    ws.ds_cnt = 1
    @inbounds ws.d[Int(src)] = zero(W)
    ws.dirty_d[1] = Int(src)
    
    ws.dh_cnt = 0 # Ensure heap dirty count is reset

    # Calculate initial bound
    log2_n1 = log2(W(n + 1))
    B = g.mean_weight * log2_n1 * W(4.0)

    # Call BMSSP with source passed in a temporary single-element buffer
    # Task local storage can be used for this too, but for one element it's fine.
    # To be zero-allocation, we can use a pre-allocated field in ws.
    # Re-purposing piv_bufs[1] for the initial call.
    ws.piv_bufs[1][1] = src
    bmsp_rec!(g, ws.d, ws.pr, ws.piv_bufs[1], 0, 1, B, 0, ws)

    return copy(ws.d), copy(ws.pr)
end

function verify_correctness()
    g = CSRGraph(4, [(1,2,2.0),(1,3,3.0),(2,4,1.0),(3,4,1.0)])
    d1, _ = ssp_duan(g, 1)
    dr, _ = dijkstra_ref(g, 1)
    return d1 ≈ dr
end

end # module
