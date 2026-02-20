module DijkstraModule

using ..CSRGraphModule: CSRGraph
using ..Common

export dijkstra_ref

mutable struct DijkstraWorkspace{T, W}
    d::Vector{W}; pr::Vector{T}
    h_vals::Vector{W}; h_idxs::Vector{T}; h_pos::Vector{T}
    dirty_d::Vector{Int}; ds_cnt::Int
end

function get_dijkstra_ws(n::Int, ::Type{T}, ::Type{W}) where {T, W}
    tls = task_local_storage(); key = (:dijkstra_v560_ws, T, W)
    ws = get!(tls, key) do
        DijkstraWorkspace{T, W}(fill(typemax(W), n), zeros(T, n), Vector{W}(undef, n), Vector{T}(undef, n), zeros(T, n), zeros(Int, n), 0)
    end::DijkstraWorkspace{T, W}
    if length(ws.d) != n
        ws = DijkstraWorkspace{T, W}(fill(typemax(W), n), zeros(T, n), Vector{W}(undef, n), Vector{T}(undef, n), zeros(T, n), zeros(Int, n), 0)
        tls[key] = ws
    end
    return ws
end

function dijkstra_ref(g::CSRGraph{T,W}, src::T) where {T,W}
    n = g.n; ws = get_dijkstra_ws(Int(n), T, W)
    
    # Selective reset
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
    
    ws.ds_cnt = 1; @inbounds ws.d[src] = zero(W); ws.dirty_d[1] = Int(src)

    fill!(ws.h_pos, zero(T)); h = Fast4AryHeap{T,W}(ws.h_vals, ws.h_idxs, ws.h_pos)
    sz = push_dec!(h, 0, src, zero(W))
    adj, wts, off = g.adjacency, g.weights, g.offset
    while sz > 0
        du, u, sz = pop_min!(h, sz); u_i = Int(u)
        @inbounds du > ws.d[u_i] && continue
        @inbounds si, ei = off[u_i], off[u_i+1]-1
        @inbounds for i in si:ei
            v, w = adj[i], wts[i]; v_i = Int(v)
            if du + w < ws.d[v_i]
                if ws.d[v_i] == typemax(W); ws.ds_cnt += 1; ws.dirty_d[ws.ds_cnt] = v_i end
                ws.d[v_i], ws.pr[v_i] = du + w, u
                sz = push_dec!(h, sz, v, ws.d[v_i])
            end
        end
    end
    return copy(ws.d), copy(ws.pr)
end

end # module

