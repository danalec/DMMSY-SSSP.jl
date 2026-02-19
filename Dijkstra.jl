"""
Dijkstra.jl — High-Performance Reference.
"""
module DijkstraModule

using ..CSRGraphModule: CSRGraph

export dijkstra_ref

struct Fast4AryHeap{T, W}
    vals::Vector{W}; idxs::Vector{T}; pos::Vector{T}
end

@inline function push_dec!(h::Fast4AryHeap{T,W}, sz::Int, n::T, d::W) where {T,W}
    @inbounds p = h.pos[Int(n)]
    if p == 0; sz += 1; i = sz
    else; i = Int(p); @inbounds if d >= h.vals[i]; return sz end
    end
    while i > 1
        par = (i - 2) >>> 2 + 1; @inbounds pv = h.vals[par]; pv <= d && break
        @inbounds pn = h.idxs[par]; h.vals[i] = pv; h.idxs[i] = pn; h.pos[Int(pn)] = i; i = par
    end
    @inbounds h.vals[i] = d; h.idxs[i] = n; h.pos[Int(n)] = i
    return sz
end

@inline function pop_min!(h::Fast4AryHeap{T,W}, sz::Int) where {T,W}
    @inbounds mv, mn = h.vals[1], h.idxs[1]; h.pos[Int(mn)] = typemax(T)
    if sz == 1 return mv, mn, 0 end
    @inbounds lv, ln = h.vals[sz], h.idxs[sz]; sz -= 1; i = 1
    while true
        c1 = (i << 2) - 2; c1 > sz && break
        mc, mcv = c1, h.vals[c1]
        c2 = c1 + 1; @inbounds if c2 <= sz && h.vals[c2] < mcv; mc, mcv = c2, h.vals[c2] end
        c3 = c1 + 2; @inbounds if c3 <= sz && h.vals[c3] < mcv; mc, mcv = c3, h.vals[c3] end
        c4 = c1 + 3; @inbounds if c4 <= sz && h.vals[c4] < mcv; mc, mcv = c4, h.vals[c4] end
        lv <= mcv && break
        @inbounds cn = h.idxs[mc]; h.vals[i] = mcv; h.idxs[i] = cn; h.pos[Int(cn)] = i; i = mc
    end
    @inbounds h.vals[i] = lv; h.idxs[i] = ln; h.pos[Int(ln)] = i
    return mv, mn, sz
end

mutable struct DijkstraWorkspace{T, W}
    d::Vector{W}; pr::Vector{T}
    h_vals::Vector{W}; h_idxs::Vector{T}; h_pos::Vector{T}
    dirty_d::Vector{Int}; ds_cnt::Int
end

function get_dijkstra_ws(n::Int, ::Type{T}, ::Type{W}) where {T, W}
    tls = task_local_storage(); key = (:dijkstra_v510_ws, T, W)
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
    if ws.ds_cnt > (n >> 2); fill!(ws.d, typemax(W))
    else; @inbounds for i in 1:ws.ds_cnt; ws.d[ws.dirty_d[i]] = typemax(W) end
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
