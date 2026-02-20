module CSRGraphModule

export CSRGraph, random_graph, reconstruct_path

struct CSRGraph{T<:Integer, W<:Real}
    n::T; m::T
    offset::Vector{T}; adjacency::Vector{T}; weights::Vector{W}
    mean_weight::W
end

function CSRGraph{T, W}(n::Int, edges::Vector{Tuple{Int, Int, W}}) where {T<:Integer, W<:Real}
    # Count edges per source vertex to determine offsets
    counts = zeros(Int, n)
    for (u, v, w) in edges
        if 1 <= u <= n && 1 <= v <= n && w >= 0
            counts[u] += 1
        end
    end
    
    m = sum(counts)
    off = Vector{T}(undef, n + 1)
    adj = Vector{T}(undef, m)
    wts = Vector{W}(undef, m)
    
    # Prefix sum to compute offsets
    off[1] = 1
    for i in 1:n
        off[i+1] = off[i] + counts[i]
    end
    
    # Reuse counts as current index pointers
    fill!(counts, 0)
    for (u, v, w) in edges
        if 1 <= u <= n && 1 <= v <= n && w >= 0
            idx = off[u] + counts[u]
            adj[idx] = T(v)
            wts[idx] = w
            counts[u] += 1
        end
    end
    
    mw = m > 0 ? sum(wts) / m : zero(W)
    return CSRGraph{T, W}(T(n), T(m), off, adj, wts, mw)
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
