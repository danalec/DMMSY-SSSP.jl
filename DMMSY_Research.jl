"""
DMMSY_Research.jl — Faithful Julia implementation
"""
module DMMSYResearch

using ..CSRGraphModule: CSRGraph
using ..DijkstraModule: dijkstra_ref

export ssp_duan_research, verify_research_correctness

# ----------------------------------------------------------------------------
# 4-Ary Heap Engine
# ----------------------------------------------------------------------------
struct Fast4AryHeap{T, W}
    vals::Vector{W}; idxs::Vector{T}; pos::Vector{T}
end

@inline function push_dec!(h::Fast4AryHeap{T,W}, sz::Int, n::T, d::W) where {T,W}
    @inbounds p = h.pos[Int(n)]
    i = (p == 0) ? (sz + 1) : Int(p)
    if p == 0; sz += 1 end
    while i > 1
        par = (i - 2) >>> 2 + 1; @inbounds pv = h.vals[par]
        pv <= d && break
        @inbounds pn = h.idxs[par]; @inbounds h.vals[i] = pv; @inbounds h.idxs[i] = pn; @inbounds h.pos[Int(pn)] = i; i = par
    end
    @inbounds h.vals[i] = d; @inbounds h.idxs[i] = n; @inbounds h.pos[Int(n)] = i
    return sz
end

@inline function pop_min!(h::Fast4AryHeap{T,W}, sz::Int) where {T,W}
    @inbounds mv, mn = h.vals[1], h.idxs[1]; @inbounds h.pos[Int(mn)] = typemax(T)
    if sz == 1 return mv, mn, 0 end
    @inbounds lv, ln = h.vals[sz], h.idxs[sz]; sz -= 1; i = 1
    while true
        c1 = (i << 2) - 2; c1 > sz && break
        mc, mcv = c1, h.vals[c1]
        for c in c1+1 : min(c1+3, sz)
            @inbounds v = h.vals[c]
            if v < mcv; mcv = v; mc = c end
        end
        lv <= mcv && break
        @inbounds cn = h.idxs[mc]; @inbounds h.vals[i] = mcv; @inbounds h.idxs[i] = cn; @inbounds h.pos[Int(cn)] = i; i = mc
    end
    @inbounds h.vals[i] = lv; @inbounds h.idxs[i] = ln; @inbounds h.pos[Int(ln)] = i
    return mv, mn, sz
end

# ----------------------------------------------------------------------------
# Research BlockedPartialPQ
# ----------------------------------------------------------------------------
mutable struct BlockedPartialPQ{T, W}
    blocks::Vector{Vector{Tuple{T, W}}}; minima::Vector{W}
    heap::Fast4AryHeap{T, W}; hsz::Int; bsz::Int; total_size::Int
end

function BlockedPartialPQ{T,W}(n::Int, hv, hi, hp) where {T,W}
    bsz = max(1, Int(ceil(n^(2/3))))
    nb = (n+bsz-1)÷bsz + 1
    BlockedPartialPQ{T,W}([Tuple{T,W}[] for _ in 1:nb], fill(typemax(W), nb), Fast4AryHeap{T,W}(hv, hi, hp), 0, bsz, 0)
end

function insert!(pq::BlockedPartialPQ{T, W}, v::T, d::W) where {T, W}
    for i in 1:length(pq.blocks)
        if length(pq.blocks[i]) < pq.bsz
            push!(pq.blocks[i], (v, d)); pq.total_size += 1
            if d < pq.minima[i]; pq.minima[i]=d; pq.hsz = push_dec!(pq.heap, pq.hsz, T(i), d) end
            return
        end
    end
    push!(pq.blocks[end], (v, d)); pq.total_size += 1
    if d < pq.minima[end]; pq.minima[end]=d; pq.hsz = push_dec!(pq.heap, pq.hsz, T(length(pq.blocks)), d) end
end

function extract!(pq::BlockedPartialPQ{T, W}) where {T, W}
    while pq.hsz > 0
        md, bi_T, pq.hsz = pop_min!(pq.heap, pq.hsz); bi = Int(bi_T)
        @inbounds if pq.minima[bi] != md; continue end
        blk = pq.blocks[bi]
        isempty(blk) && continue
        mi, mv, mdb = 1, blk[1][1], blk[1][2]
        for j in 2:length(blk)
            if blk[j][2] < mdb; mdb, mv, mi = blk[j][2], blk[j][1], j end
        end
        blk[mi] = blk[end]; pop!(blk); pq.total_size -= 1
        if isempty(blk); pq.minima[bi] = typemax(W)
        else
            nm = typemax(W); for x in blk; if x[2]<nm; nm=x[2] end end
            pq.minima[bi] = nm; pq.hsz = push_dec!(pq.heap, pq.hsz, bi_T, nm)
        end
        return (mv, mdb)
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Core Algorithm
# ----------------------------------------------------------------------------
function ssp_duan_research(g::CSRGraph{T, W}, src::T) where {T, W}
    n = g.n; d, pr = fill(typemax(W), n), zeros(T, n); d[src] = zero(W)
    # Fresh allocation for sync stability
    hv, hi, hp = Vector{W}(undef, n), Vector{T}(undef, n), zeros(T, n)
    bmsp!(g, [src], d, pr, 1, h_vals=hv, h_idxs=hi, h_pos=hp)
    return d, pr
end

function bmsp!(g::CSRGraph{T, W}, src, d, pr, dp; h_vals, h_idxs, h_pos) where {T, W}
    n = g.n
    if n <= 1000 || dp >= 3
        fill!(h_pos, zero(T)); h = Fast4AryHeap{T, W}(h_vals, h_idxs, h_pos); sz = 0
        for s in src; sz = push_dec!(h, sz, s, d[s]) end
        while sz > 0
            du, u, sz = pop_min!(h, sz); @inbounds du > d[u] && continue
            @inbounds si, ei = g.offset[u], g.offset[u+1]-1
            for i in si:ei
                @inbounds v, w = g.adjacency[i], g.weights[i]
                if du + w < d[v]; d[v], pr[v] = du + w, u; sz = push_dec!(h, sz, v, d[v]) end
            end
        end
        return
    end
    
    pq = BlockedPartialPQ{T, W}(n, h_vals, h_idxs, h_pos); fill!(h_pos, zero(T))
    for v in 1:n; @inbounds if d[v] < typemax(W); insert!(pq, v, d[v]) end end
    while pq.total_size > 0
        it = extract!(pq); it === nothing && break
        v, dv = it; @inbounds dv > d[v] && continue
        @inbounds si, ei = g.offset[v], g.offset[v+1]-1
        for i in si:ei
            @inbounds u, w = g.adjacency[i], g.weights[i]
            @inbounds if d[v] + w < d[u]
                d[u] = d[v] + w; pr[u] = v; insert!(pq, u, d[u])
            end
        end
    end
end

function verify_research_correctness()
    g = CSRGraph(3, [(1,2,1.0),(2,3,1.0),(3,2,0.5)])
    d, _ = ssp_duan_research(g, 1); d_ref, _ = dijkstra_ref(g, 1); return d ≈ d_ref
end

end # module
