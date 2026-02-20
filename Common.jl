
module Common
export Fast4AryHeap, push_dec!, pop_min!

# ----------------------------------------------------------------------------
# Optimized 4-Ary Heap used across multiple SSSP implementations
# ----------------------------------------------------------------------------
struct Fast4AryHeap{T<:Integer, W<:Real}
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

end # module
