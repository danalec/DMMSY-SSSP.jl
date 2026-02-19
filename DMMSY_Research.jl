"""
DMMSY_Research.jl — Simplified DMMSY Implementation
Based on the key ideas of "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths"
"""
module DMMSYResearch

using ..CSRGraphModule: CSRGraph
using ..DijkstraModule: dijkstra_ref

export ssp_duan_research, verify_research_correctness

# ----------------------------------------------------------------------------
# Algorithm Parameters (using log2 for efficiency)
# ----------------------------------------------------------------------------
function get_params(n::Int)
    k = max(4, Int(floor(log2(n)^(1/3))))
    t = max(2, Int(floor(log2(n)^(2/3))))
    return k, t
end

# ----------------------------------------------------------------------------
# 4-Ary Heap for Base Case
# ----------------------------------------------------------------------------
struct Fast4AryHeap{T, W}
    vals::Vector{W}
    idxs::Vector{T}
    pos::Vector{T}
end

@inline function push_dec!(h::Fast4AryHeap{T,W}, sz::Int, n::T, d::W) where {T,W}
    @inbounds p = h.pos[Int(n)]
    if p == 0
        sz += 1
        i = sz
    else
        i = Int(p)
        @inbounds if d >= h.vals[i]; return sz end
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
    return sz
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
        c1 > sz && break
        mc, mcv = c1, h.vals[c1]
        c2 = c1 + 1
        @inbounds if c2 <= sz && h.vals[c2] < mcv; mc, mcv = c2, h.vals[c2] end
        c3 = c1 + 2
        @inbounds if c3 <= sz && h.vals[c3] < mcv; mc, mcv = c3, h.vals[c3] end
        c4 = c1 + 3
        @inbounds if c4 <= sz && h.vals[c4] < mcv; mc, mcv = c4, h.vals[c4] end
        lv <= mcv && break
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
# Optimized Workspace with pre-allocated buffers
# ----------------------------------------------------------------------------
mutable struct ResearchWorkspace{T, W}
    d::Vector{W}
    pr::Vector{T}
    h_vals::Vector{W}
    h_idxs::Vector{T}
    h_pos::Vector{T}
    dirty_d::Vector{Int}
    ds_cnt::Int
    piv_buf::Vector{T}
    active_buf::Vector{T}
end

function get_research_ws(n::Int, ::Type{T}, ::Type{W}) where {T, W}
    tls = task_local_storage()
    key = (:dmmsy_research_ws, T, W)
    ws = get!(tls, key) do
        ResearchWorkspace{T, W}(
            fill(typemax(W), n),
            zeros(T, n),
            Vector{W}(undef, n),
            Vector{T}(undef, n),
            zeros(T, n),
            zeros(Int, n),
            0,
            Vector{T}(undef, n),
            Vector{T}(undef, n)
        )
    end::ResearchWorkspace{T, W}
    if length(ws.d) != n
        ws = ResearchWorkspace{T, W}(
            fill(typemax(W), n),
            zeros(T, n),
            Vector{W}(undef, n),
            Vector{T}(undef, n),
            zeros(T, n),
            zeros(Int, n),
            0,
            Vector{T}(undef, n),
            Vector{T}(undef, n)
        )
        tls[key] = ws
    end
    return ws
end

# ----------------------------------------------------------------------------
# Optimized BMSSP Implementation
# ----------------------------------------------------------------------------
function bmsp!(g::CSRGraph{T, W}, d::Vector{W}, pr::Vector{T}, src::Vector{T}, B::W, dp::Int, ws::ResearchWorkspace{T, W}) where {T, W}
    n = g.n
    k, t = get_params(n)
    off, adj, wts = g.offset, g.adjacency, g.weights

    # Base case: dp >= t or small src (no bound filtering for speed)
    if dp >= t || length(src) <= k
        fill!(ws.h_pos, zero(T))
        h = Fast4AryHeap{T, W}(ws.h_vals, ws.h_idxs, ws.h_pos)
        sz = 0

        for s in src
            sz = push_dec!(h, sz, s, d[Int(s)])
        end

        while sz > 0
            du, u, sz = pop_min!(h, sz)
            u_i = Int(u)
            @inbounds if du > d[u_i]; continue end

            @inbounds si, ei = off[u_i], off[u_i+1]-1
            for i in si:ei
                @inbounds v, w = adj[i], wts[i]
                v_i = Int(v)
                nd = du + w
                if nd <= d[v_i]
                    if d[v_i] == typemax(W)
                        ws.ds_cnt += 1
                        ws.dirty_d[ws.ds_cnt] = v_i
                    end
                    d[v_i] = nd
                    pr[v_i] = u
                    sz = push_dec!(h, sz, v, nd)
                end
            end
        end

        return
    end

    # Select k pivots evenly from src using pre-allocated buffer (zero allocation!)
    np = min(length(src), k)
    piv_buf = ws.piv_buf
    step = max(1, length(src) ÷ np)
    np = 0
    for i in 1:step:length(src)
        np += 1
        piv_buf[np] = src[i]
        if np >= k; break end
    end
    pivots = collect(piv_buf[1:np])

    # Recursive call on pivots with tighter bound
    bmsp!(g, d, pr, pivots, B * W(0.5), dp + 1, ws)

    # Main loop - process remaining src vertices with bound B
    fill!(ws.h_pos, zero(T))
    h = Fast4AryHeap{T, W}(ws.h_vals, ws.h_idxs, ws.h_pos)
    sz = 0

    for s in src
        s_i = Int(s)
        @inbounds if d[s_i] < B
            sz = push_dec!(h, sz, s, d[s_i])
        end
    end

    while sz > 0
        du, u, sz = pop_min!(h, sz)
        u_i = Int(u)
        @inbounds if du > d[u_i]; continue end

        @inbounds si, ei = off[u_i], off[u_i+1]-1
        for i in si:ei
            @inbounds v, w = adj[i], wts[i]
            v_i = Int(v)
            nd = du + w
            if nd <= d[v_i]
                if d[v_i] == typemax(W)
                    ws.ds_cnt += 1
                    ws.dirty_d[ws.ds_cnt] = v_i
                end
                d[v_i] = nd
                pr[v_i] = u
                @inbounds if nd < B
                    sz = push_dec!(h, sz, v, nd)
                end
            end
        end
    end
end

# ----------------------------------------------------------------------------
# Main Interface Function
# ----------------------------------------------------------------------------
function ssp_duan_research(g::CSRGraph{T, W}, src::T) where {T, W}
    n = g.n
    if n == 0
        return W[], T[]
    end

    ws = get_research_ws(n, T, W)

    # Initialize distances - use selective reset
    if ws.ds_cnt > (n >> 2)
        fill!(ws.d, typemax(W))
    else
        @inbounds for i in 1:ws.ds_cnt
            ws.d[ws.dirty_d[i]] = typemax(W)
        end
    end

    ws.ds_cnt = 1
    @inbounds ws.d[Int(src)] = zero(W)
    ws.dirty_d[1] = Int(src)

    # Calculate initial bound using log2 and larger multiplier for better performance
    B = g.mean_weight * log2(n + 1) * W(4.0)

    # Call BMSSP with active buffer for source
    ws.active_buf[1] = src
    bmsp!(g, ws.d, ws.pr, collect(ws.active_buf[1:1]), B, 0, ws)

    return copy(ws.d), copy(ws.pr)
end

function verify_research_correctness()
    g = CSRGraph(3, [(1,2,1.0),(2,3,1.0),(3,2,0.5)])
    d, _ = ssp_duan_research(g, 1)
    d_ref, _ = dijkstra_ref(g, 1)
    return d ≈ d_ref
end

end # module
