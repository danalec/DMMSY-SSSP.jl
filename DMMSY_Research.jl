module DMMSYResearch

using ..CSRGraphModule: CSRGraph
using ..DijkstraModule: dijkstra_ref
using ..Common

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
    key = (:dmmsy_research_v560_ws, T, W)
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
function bmsp!(g::CSRGraph{T, W}, d::Vector{W}, pr::Vector{T}, src::AbstractVector{T}, B::W, dp::Int, ws::ResearchWorkspace{T, W}) where {T, W}
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
    step = max(1, length(src) / np)
    curr_np = 0
    for i in 1:step:length(src)
        idx = Int(floor(i))
        if idx > length(src); break end
        curr_np += 1
        piv_buf[curr_np] = src[idx]
        if curr_np >= k; break end
    end

    # Recursive call on pivots with tighter bound
    bmsp!(g, d, pr, view(piv_buf, 1:curr_np), B * W(0.5), dp + 1, ws)

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

    # Calculate initial bound
    B = g.mean_weight * log2(n + 1) * W(4.0)

    # Call BMSSP with active buffer for source
    ws.active_buf[1] = src
    bmsp!(g, ws.d, ws.pr, view(ws.active_buf, 1:1), B, 0, ws)

    return copy(ws.d), copy(ws.pr)
end

end # module

