module DMMSYResearch

using ..CSRGraphModule: CSRGraph
using ..DijkstraModule: dijkstra_ref
using ..Common: get_params, Edge, HeapNode, push_dec!, pop_min!, Fast4AryHeap

export ssp_duan_research, verify_research_correctness

# ----------------------------------------------------------------------------
# Optimized Workspace with pre-allocated buffers
# ----------------------------------------------------------------------------
mutable struct ResearchWorkspace{T, W}
    d::Vector{W}
    pr::Vector{T}
    h_nodes::Vector{HeapNode{T, W}}
    h_pos::Vector{T}
    dirty_h::Vector{T}
    dirty_d::Vector{Int}
    ds_cnt::Int
    piv_bufs::Vector{Vector{T}}
    active_buf::Vector{T}
end

function get_research_ws(n::Int, ::Type{T}, ::Type{W}) where {T, W}
    tls = task_local_storage()
    key = (:dmmsy_research_v560_ws, T, W)
    ws = get!(tls, key) do
        ResearchWorkspace{T, W}(
            fill(typemax(W), n),
            zeros(T, n),
            Vector{HeapNode{T, W}}(undef, n),
            zeros(T, n),
            zeros(T, n),
            zeros(Int, n),
            0,
            [Vector{T}(undef, max(4, Int(floor(log2(n)^(1/3))))) for _ in 1:max(2, Int(floor(log2(n)^(2/3)))) + 1],
            Vector{T}(undef, n)
        )
    end::ResearchWorkspace{T, W}
    if length(ws.d) != n
        ws = ResearchWorkspace{T, W}(
            fill(typemax(W), n),
            zeros(T, n),
            Vector{HeapNode{T, W}}(undef, n),
            zeros(T, n),
            zeros(T, n),
            zeros(Int, n),
            0,
            [Vector{T}(undef, max(4, Int(floor(log2(n)^(1/3))))) for _ in 1:max(2, Int(floor(log2(n)^(2/3)))) + 1],
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
    off, edges = g.offset, g.edges

    # Base case: dp >= t or small src (no bound filtering for speed)
    if dp >= t || length(src) <= k
        fill!(ws.h_pos, zero(T))
        h = Fast4AryHeap{T, W}(ws.h_nodes, ws.h_pos, ws.dirty_h)
        sz, dcnt = 0, 0

        for s in src
            sz, dcnt = push_dec!(h, sz, dcnt, s, d[Int(s)])
        end

        while sz > 0
            du, u, sz = pop_min!(h, sz)
            u_i = Int(u)
            @inbounds if du > d[u_i]; continue end

            @inbounds si, ei = off[u_i], off[u_i+1]-1
            for i in si:ei
                @inbounds e = edges[i]
                @inbounds v, w = e.v, e.w
                v_i = Int(v)
                nd = du + w
                if nd <= d[v_i]
                    if d[v_i] == typemax(W)
                        ws.ds_cnt += 1
                        ws.dirty_d[ws.ds_cnt] = v_i
                    end
                    d[v_i] = nd
                    pr[v_i] = u
                    sz, dcnt = push_dec!(h, sz, dcnt, v, nd)
                end
            end
        end

        return
    end

    # Select k pivots evenly from src using pre-allocated buffer (zero allocation!)
    np = min(length(src), k)
    pivots = ws.piv_bufs[dp + 2]
    step = max(1, length(src) ÷ np)
    curr_np = 0
    bound = min(length(src), step * k)
    for i in 1:step:bound
        curr_np += 1
        @inbounds pivots[curr_np] = src[i]
    end

    # Recursive call on pivots with tighter bound
    bmsp!(g, d, pr, view(pivots, 1:curr_np), B * W(0.5), dp + 1, ws)

    # Main loop - process remaining src vertices with bound B
    fill!(ws.h_pos, zero(T))
    h = Fast4AryHeap{T, W}(ws.h_nodes, ws.h_pos, ws.dirty_h)
    sz, dcnt = 0, 0

    for s in src
        s_i = Int(s)
        @inbounds if d[s_i] < B
            sz, dcnt = push_dec!(h, sz, dcnt, s, d[s_i])
        end
    end

    while sz > 0
        du, u, sz = pop_min!(h, sz)
        u_i = Int(u)
        @inbounds if du > d[u_i]; continue end

        @inbounds si, ei = off[u_i], off[u_i+1]-1
        for i in si:ei
            @inbounds e = edges[i]
            @inbounds v, w = e.v, e.w
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
                    sz, dcnt = push_dec!(h, sz, dcnt, v, nd)
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

