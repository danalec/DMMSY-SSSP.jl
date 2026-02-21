
module Common
export Fast4AryHeap, push_dec!, pop_min!, Edge, get_params, HeapNode

# ----------------------------------------------------------------------------
# Interleaved Edge structure for better cache locality
# ----------------------------------------------------------------------------
struct Edge{T<:Integer, W<:Real}
    v::T
    w::W
end

# ----------------------------------------------------------------------------
# Shared algorithm parameters
# ----------------------------------------------------------------------------
@inline function get_params(n::Integer)
    k = max(4, Int(floor(log2(n)^(1/3))))
    t = max(2, Int(floor(log2(n)^(2/3))))
    return k, t
end

# ----------------------------------------------------------------------------
# Optimized 4-Ary Heap used across multiple SSSP implementations
# ----------------------------------------------------------------------------
struct HeapNode{T<:Integer, W<:Real}
    v::W
    i::T
end

struct Fast4AryHeap{T<:Integer, W<:Real}
    nodes::Vector{HeapNode{T, W}}
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
        @inbounds if d >= h.nodes[i].v; return sz, dcnt end
    end
    while i > 1
        par = (i - 2) >>> 2 + 1
        @inbounds node = h.nodes[par]
        node.v <= d && break
        @inbounds h.nodes[i] = node
        @inbounds h.pos[Int(node.i)] = i
        i = par
    end
    @inbounds h.nodes[i] = HeapNode{T, W}(d, n)
    @inbounds h.pos[Int(n)] = i
    return sz, dcnt
end

@inline function pop_min!(h::Fast4AryHeap{T,W}, sz::Int) where {T,W}
    @inbounds mn_node = h.nodes[1]
    mv, mn = mn_node.v, mn_node.i
    @inbounds h.pos[Int(mn)] = typemax(T)
    if sz == 1; return mv, mn, 0 end
    @inbounds last_node = h.nodes[sz]
    lv, ln = last_node.v, last_node.i
    sz -= 1
    i = 1
    while true
        c1 = (i << 2) - 2
        c1 > sz && break
        @inbounds mc_node = h.nodes[c1]
        mc, mcv = c1, mc_node.v
        
        @inbounds if c1 + 1 <= sz
            c2_node = h.nodes[c1 + 1]
            if c2_node.v < mcv; mc, mcv, mc_node = c1 + 1, c2_node.v, c2_node end
        end
        @inbounds if c1 + 2 <= sz
            c3_node = h.nodes[c1 + 2]
            if c3_node.v < mcv; mc, mcv, mc_node = c1 + 2, c3_node.v, c3_node end
        end
        @inbounds if c1 + 3 <= sz
            c4_node = h.nodes[c1 + 3]
            if c4_node.v < mcv; mc, mcv, mc_node = c1 + 3, c4_node.v, c4_node end
        end
        
        lv <= mcv && break
        
        @inbounds h.nodes[i] = mc_node
        @inbounds h.pos[Int(mc_node.i)] = i
        i = mc
    end
    @inbounds h.nodes[i] = last_node
    @inbounds h.pos[Int(ln)] = i
    return mv, mn, sz
end

end # module
