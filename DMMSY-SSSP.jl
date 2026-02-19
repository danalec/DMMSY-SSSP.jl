"""
DMMSY-SSSP.jl — State-of-the-art SSSP.
"""
module DMMSYSSSP

using Printf
using ..CSRGraphModule: CSRGraph, random_graph
using ..DijkstraModule: dijkstra_ref

export ssp_duan, verify_correctness

# ----------------------------------------------------------------------------
# Ultra-Tight 4-Ary Heap
# ----------------------------------------------------------------------------
struct Fast4AryHeap{T<:Integer, W<:Real}
    vals::Vector{W}; idxs::Vector{T}; pos::Vector{T}
end

@inline function push_dec!(h::Fast4AryHeap{T,W}, sz::Int, n::T, d::W) where {T,W}
    @inbounds p = h.pos[Int(n)]
    if p == 0
        sz += 1; i = sz
    else
        i = Int(p)
        @inbounds if d >= h.vals[i]; return sz end
    end
    while i > 1
        par = (i - 2) >>> 2 + 1; @inbounds pv = h.vals[par]
        pv <= d && break
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

# ----------------------------------------------------------------------------
# Workspace and Recursive Structure
# ----------------------------------------------------------------------------
struct Layer{T, W}
    h_vals::Vector{W}; h_idxs::Vector{T}; h_pos::Vector{T}
end

mutable struct WorkSpace{T<:Integer, W<:Real}
    d::Vector{W}; pr::Vector{T}
    data::Vector{W}; bm::Vector{UInt64}
    bh_vals::Vector{W}; bh_idxs::Vector{T}; bh_pos::Vector{T}
    layers::Vector{Layer{T, W}}
    dirty_data::Vector{Int}; dd_cnt::Int
    dirty_d::Vector{Int}; ds_cnt::Int
    active_buf::Vector{T}
    piv_buf::Vector{T}
end

function WorkSpace{T, W}(n::Int, m::Int) where {T, W}
    bb = 6; bsz = 1 << bb; nb = (n+bsz-1)>>bb; nm = (nb+63)>>6
    layers = [Layer{T, W}(Vector{W}(undef, n), Vector{T}(undef, n), zeros(T, n)) for _ in 1:4]
    WorkSpace{T, W}(fill(typemax(W), n), zeros(T, n),
                    fill(typemax(W), nb*(bsz+1)), zeros(UInt64, nm),
                    Vector{W}(undef, nb), Vector{T}(undef, nb), zeros(T, nb),
                    layers, Vector{Int}(undef, nb*(bsz+1)), 0, Vector{Int}(undef, n), 0,
                    Vector{T}(undef, n), Vector{T}(undef, n))
end

function get_workspace(n::Int, m::Int, ::Type{T}, ::Type{W}) where {T, W}
    tls = task_local_storage(); key = (:dmmsy_v510_ws, T, W)
    ws = get!(tls, key) do; WorkSpace{T, W}(n, m) end::WorkSpace{T, W}
    if length(ws.d) < n; ws = WorkSpace{T, W}(n, m); tls[key] = ws end
    return ws
end

function ssp_duan(g::CSRGraph{T, W}, src::T) where {T, W}
    n = g.n; n == 0 && return W[], T[]
    ws = get_workspace(n, g.m, T, W)
    
    # Selective d-reset
    if ws.ds_cnt > (n >> 2); fill!(ws.d, typemax(W))
    else; @inbounds for i in 1:ws.ds_cnt; ws.d[ws.dirty_d[i]] = typemax(W) end
    end
    ws.ds_cnt = 0; @inbounds ws.d[src] = zero(W); ws.ds_cnt += 1; ws.dirty_d[1] = Int(src)
    
    th = g.mean_weight * log2(n+1) * 3.0
    bmsp_rec!(g, T[src], ws.d, ws.pr, W(th), 1, ws, 6, 64)
    
    # Return copies as mandated for safety
    return copy(ws.d), copy(ws.pr)
end

@inline function record_data!(ws, idx)
    ws.dd_cnt += 1; @inbounds ws.dirty_data[ws.dd_cnt] = idx
end

@inline function record_dist!(ws, idx)
    ws.ds_cnt += 1; @inbounds ws.dirty_d[ws.ds_cnt] = idx
end

function bmsp_rec!(g::CSRGraph{T,W}, src::Vector{T}, d::Vector{W}, pr::Vector{T}, th::W, dp::Int, ws::WorkSpace{T,W}, bb, bsz) where {T,W}
    if dp >= 3 || length(src) < 32
        h = Fast4AryHeap(ws.layers[dp].h_vals, ws.layers[dp].h_idxs, ws.layers[dp].h_pos)
        fill!(h.pos, zero(T)); base_case!(g, d, pr, src, h, ws); return
    end
    
    # Picking pivots from src + a small sample of active frontier
    piv = T[]; push!(piv, src[1])
    np = min(length(src), 15); for i in 1:max(1, length(src)÷np):length(src); push!(piv, src[i]); length(piv)>=np && break end
    
    bmsp_rec!(g, piv, d, pr, th*W(0.4), dp+1, ws, bb, bsz)
    propagate!(ws, g, src, d, pr, bb, bsz)
end

function propagate!(ws::WorkSpace{T,W}, g::CSRGraph{T,W}, sources::Vector{T}, d::Vector{W}, pr::Vector{T}, bb, bsz) where {T,W}
    data, bm, bh_vals, bh_idxs, bh_pos = ws.data, ws.bm, ws.bh_vals, ws.bh_idxs, ws.bh_pos
    if ws.dd_cnt > (length(data) >> 3); fill!(data, typemax(W))
    else; @inbounds for i in 1:ws.dd_cnt; data[ws.dirty_data[i]] = typemax(W) end
    end
    ws.dd_cnt = 0; fill!(bm, 0); fill!(bh_pos, zero(T)); hsz = 0; nb = length(bh_vals)
    
    @inbounds for s in sources
        s_i = Int(s); idx = ((s_i-1)>>bb)*(bsz+1) + ((s_i-1)&(bsz-1)) + 2
        ds = d[s_i]; data[idx] = ds; record_data!(ws, idx)
        b = ((s_i-1)>>bb)+1; mi = (b-1)*(bsz+1)+1
        if ds < data[mi]
            data[mi] = ds; record_data!(ws, mi)
            bm[(b-1)>>6+1] |= (UInt64(1)<<((b-1)&63))
            if nb > 256; hsz = push_dec!(Fast4AryHeap(bh_vals, bh_idxs, bh_pos), hsz, T(b), ds) end
        end
    end
    
    off, adj, wts = g.offset, g.adjacency, g.weights; n = g.n
    while true
        bb_id, bd = 0, typemax(W)
        if nb <= 256
            @inbounds for m in 1:length(bm)
                mv = bm[m]; mv == 0 && continue
                while mv != 0
                    bi = trailing_zeros(mv); b = ((m - 1) << 6) + bi + 1
                    dist = data[(b - 1) * (bsz + 1) + 1]
                    if dist < bd; bd, bb_id = dist, b end
                    mv &= (mv - 1)
                end
            end
        else
            h = Fast4AryHeap(bh_vals, bh_idxs, bh_pos)
            while hsz > 0
                val, b_id_T, hsz = pop_min!(h, hsz); b = Int(b_id_T)
                @inbounds if data[(b-1)*(bsz+1)+1] == val; bb_id, bd = b, val; break end
            end
        end
        bb_id == 0 && break
        
        sv = ((bb_id-1)<<bb)+1; ev = min(sv+bsz-1, n); u, ud = zero(T), typemax(W); ofst = (bb_id-1)*(bsz+1)
        @inbounds for i in 1:(ev-sv+1)
            dv = data[ofst+1+i]
            if dv < ud; ud, u = dv, T(sv+i-1) end
        end
        
        if u == 0; @inbounds bm[(bb_id-1)>>6+1] &= ~(UInt64(1)<<((bb_id-1)&63)); continue end
        
        u_i = Int(u); iu = ((u_i-1)>>bb)*(bsz+1) + ((u_i-1)&(bsz-1)) + 2
        @inbounds data[iu] = typemax(W)
        nm = typemax(W); @inbounds for i in 1:(ev-sv+1); nm = min(nm, data[ofst+1+i]) end
        @inbounds data[ofst+1] = nm
        
        if nm == typemax(W); @inbounds bm[(bb_id-1)>>6+1] &= ~(UInt64(1)<<((bb_id-1)&63))
        elseif nb > 256; hsz = push_dec!(Fast4AryHeap(bh_vals, bh_idxs, bh_pos), hsz, T(bb_id), nm) end
        
        @inbounds si, ei = off[u_i], off[u_i+1]-1
        @inbounds for j in si:ei
            v, w = adj[j], wts[j]; nd = ud + w
            v_i = Int(v)
            if nd < d[v_i]
                if d[v_i] == typemax(W); record_dist!(ws, v_i) end
                d[v_i], pr[v_i] = nd, u
                iv = ((v_i-1)>>bb)*(bsz+1) + ((v_i-1)&(bsz-1)) + 2
                if nd < data[iv]
                    data[iv] = nd; record_data!(ws, iv)
                    bv = ((v_i-1)>>bb)+1; mvi = (bv-1)*(bsz+1)+1
                    if nd < data[mvi]
                        data[mvi] = nd; record_data!(ws, mvi)
                        bm[(bv-1)>>6+1] |= (UInt64(1)<<((bv-1)&63))
                        if nb > 256; hsz = push_dec!(Fast4AryHeap(bh_vals, bh_idxs, bh_pos), hsz, T(bv), nd) end
                    end
                end
            end
        end
    end
end

function base_case!(g::CSRGraph{T,W}, d::Vector{W}, pr::Vector{T}, src::Vector{T}, h::Fast4AryHeap{T,W}, ws::WorkSpace) where {T,W}
    sz = 0; @inbounds for s in src; sz = push_dec!(h, sz, s, d[Int(s)]) end
    adj, wts, off = g.adjacency, g.weights, g.offset
    while sz > 0
        du, u, sz = pop_min!(h, sz); u_i = Int(u)
        @inbounds du > d[u_i] && continue
        @inbounds si, ei = off[u_i], off[u_i+1]-1
        @inbounds for i in si:ei
            v, w = adj[i], wts[i]; v_i = Int(v)
            if du + w < d[v_i]
                if d[v_i] == typemax(W); record_dist!(ws, v_i) end
                d[v_i], pr[v_i] = du + w, u
                sz = push_dec!(h, sz, v, d[v_i])
            end
        end
    end
end

function verify_correctness()
    g = CSRGraph(4, [(1,2,2.0),(1,3,3.0),(2,4,1.0),(3,4,1.0)])
    d1, _ = ssp_duan(g, 1); dr, _ = dijkstra_ref(g, 1); return d1 ≈ dr
end

end # module
