module CSRGraphModule

using ..Common: Edge

export CSRGraph, random_graph, reconstruct_path

struct CSRGraph{T<:Integer, W<:Real}
    n::T; m::T
    offset::Vector{T}
    edges::Vector{Edge{T, W}}
    mean_weight::W
end

function CSRGraph{T, W}(n::Int, edges_in::Vector{Tuple{Int, Int, W}}) where {T<:Integer, W<:Real}
    # Count edges per source vertex to determine offsets
    counts = zeros(Int, n)
    for (u, v, w) in edges_in
        if 1 <= u <= n && 1 <= v <= n && w >= 0
            counts[u] += 1
        end
    end
    
    m = sum(counts)
    off = Vector{T}(undef, n + 1)
    e_out = Vector{Edge{T, W}}(undef, m)
    
    # Prefix sum to compute offsets
    off[1] = 1
    for i in 1:n
        off[i+1] = off[i] + counts[i]
    end
    
    # Reuse counts as current index pointers
    fill!(counts, 0)
    for (u, v, w) in edges_in
        if 1 <= u <= n && 1 <= v <= n && w >= 0
            idx = off[u] + counts[u]
            e_out[idx] = Edge{T, W}(T(v), w)
            counts[u] += 1
        end
    end
    
    sum_w = zero(W)
    for i in 1:m
        @inbounds sum_w += e_out[i].w
    end
    mw = m > 0 ? sum_w / m : zero(W)
    
    return CSRGraph{T, W}(T(n), T(m), off, e_out, mw)
end

CSRGraph(n::Int, edges::Vector{Tuple{Int, Int, W}}) where W = CSRGraph{Int, W}(n, edges)

function random_graph(n::Int, m::Int, max_w::W) where W
    e = Vector{Tuple{Int, Int, W}}(undef, m)
    for i in 1:m
        e[i] = (rand(1:n), rand(1:n), rand() * max_w)
    end
    return CSRGraph{Int, W}(n, e)
end

function reconstruct_path(pred::Vector{Int}, source::Int, target::Int)
    (pred[target] == 0 && target != source) && return Int[]
    path = Int[target]; curr = target
    while curr != source && curr != 0
        curr = pred[curr]
        curr == 0 && return Int[]
        pushfirst!(path, curr)
    end
    return path
end

end
