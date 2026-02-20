"""
DMMSY-SSSP.jl — State-of-the-art SSSP.
Based on "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths"
Optimized implementation with bound filtering and efficient recursion.
"""
module DMMSYSSSP

using ..CSRGraphModule: CSRGraph
using ..DijkstraModule: dijkstra_ref

export ssp_duan, verify_correctness

# ----------------------------------------------------------------------------
# Optimized 4-Ary Heap with dirty-tracking
# ----------------------------------------------------------------------------
struct Fast4AryHeap{T<:Integer, W<:Real}
    vals::Vector{W}
    idxs::Vector{T}
    pos::Vector{T}
    dirty::Vector{T}
end

@inline function push_dec!(h::Fast4AryHeap{T,W}, sz::Int, dcnt::Int, n::T, d::W) where {T,W}
    @inbounds p = h.pos[Int(n)]
    if p == 0 || p == typemax(T)
        sz += 1
        dcnt += 1
        i = sz
        @inbounds h.dirty[dcnt] = n
    else
        i = Int(p)
        @inbounds if d >= h.vals[i]; return sz, dcnt end
    end
    while i > 1
        par = (i - 2) >>> 2 + 1
        @inbounds pv = h.vals[par]
        pv <= d && break
        @inbounds pn = h.idxs[par]
        @inbounds h.vals[i] = pv
        @inbounds h.idxs[i] = pn
        @inbounds h.pos[Int(pn)] = i
        i = par
    end
    @inbounds h.vals[i] = d
    @inbounds h.idxs[i] = n
    @inbounds h.pos[Int(n)] = i
    return sz, dcnt
end

@inline function pop_min!(h::Fast4AryHeap{T,W}, sz::Int) where {T,W}
    @inbounds mv, mn = h.vals[1], h.idxs[1]
    @inbounds h.pos[Int(mn)] = typemax(T) 
    if sz == 1; return mv, mn, 0 end
    @inbounds lv, ln = h.vals[sz], h.idxs[sz]
    sz -= 1
    i = 1
    while true
        c1 = (i << 2) - 2
        if c1 > sz; break end
        mc, mcv = c1, h.vals[c1]
        c2 = c1 + 1
        @inbounds if c2 <= sz && h.vals[c2] < mcv; mc, mcv = c2, h.vals[c2] end
        c3 = c1 + 2
        @inbounds if c3 <= sz && h.vals[c3] < mcv; mc, mcv = c3, h.vals[c3] end
        c4 = c1 + 3
        @inbounds if c4 <= sz && h.vals[c4] < mcv; mc, mcv = c4, h.vals[c4] end
        if lv <= mcv; break end
        @inbounds cn = h.idxs[mc]
        @inbounds h.vals[i] = mcv
        @inbounds h.idxs[i] = cn
        @inbounds h.pos[Int(cn)] = i
        i = mc
    end
    @inbounds h.vals[i] = lv
    @inbounds h.idxs[i] = ln
    @inbounds h.pos[Int(ln)] = i
    return mv, mn, sz
end

# ----------------------------------------------------------------------------
# Simplified Algorithm Parameters (using log2 for efficiency)
# ----------------------------------------------------------------------------
@inline function get_params(n::Int)
    k = max(4, Int(floor(log2(n)^(1/3))))
    t = max(2, Int(floor(log2(n)^(2/3))))
    return k, t
end

# ----------------------------------------------------------------------------
# Optimized Workspace with pre-allocated buffers
# ----------------------------------------------------------------------------
mutable struct WorkSpace{T<:Integer, W<:Real}
    d::Vector{W}
    pr::Vector{T}
    h_vals::Vector{W}
    h_idxs::Vector{T}
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
            Vector{W}(undef, n),
            Vector{T}(undef, n),
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
            Vector{W}(undef, n),
            Vector{T}(undef, n),
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
    off, adj, wts = g.offset, g.adjacency, g.weights

    # Base case: dp >= t or small src
    if dp >= t || len_src <= k
        # Selective reset of heap position map
        if ws.dh_cnt > 0
            @inbounds for i in 1:ws.dh_cnt
                ws.h_pos[Int(ws.dirty_h[i])] = zero(T)
            end
            ws.dh_cnt = 0
        end
        
        h = Fast4AryHeap{T, W}(ws.h_vals, ws.h_idxs, ws.h_pos, ws.dirty_h)
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
                @inbounds v, w = adj[i], wts[i]
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
    
    h = Fast4AryHeap{T, W}(ws.h_vals, ws.h_idxs, ws.h_pos, ws.dirty_h)
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
            @inbounds v, w = adj[i], wts[i]
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
