module CSRGraphModule
export CSRGraph, random_graph, reconstruct_path
struct CSRGraph{T<:Integer, W<:Real}
    n::T; m::T
    offset::Vector{T}; adjacency::Vector{T}; weights::Vector{W}
    mean_weight::W
end
function CSRGraph{T, W}(n::Int, edges::Vector{Tuple{Int, Int, W}}) where {T<:Integer, W<:Real}
    buckets = [Tuple{Int, W}[] for _ in 1:n]
    for (u, v, w) in edges
        if 1 <= u <= n && 1 <= v <= n && w >= 0; push!(buckets[u], (v, w)) end
    end
    m = sum(length(b) for b in buckets); off = zeros(T, n + 1); adj = zeros(T, m); wts = zeros(W, m)
    idx = 1
    for u in 1:n
        off[u] = idx
        for (v, w) in buckets[u]; adj[idx] = v; wts[idx] = w; idx += 1 end
    end
    off[n+1] = m + 1; mw = m > 0 ? sum(wts) / m : zero(W)
    return CSRGraph{T, W}(T(n), T(m), off, adj, wts, mw)
end
CSRGraph(n::Int, edges::Vector{Tuple{Int, Int, W}}) where W = CSRGraph{Int, W}(n, edges)

function random_graph(n::Int, m::Int, max_w::W) where W
    e = Tuple{Int, Int, W}[]
    for _ in 1:m; push!(e, (rand(1:n), rand(1:n), rand()*max_w)) end
    return CSRGraph{Int, W}(n, e)
end
function reconstruct_path(pred::Vector{Int}, source::Int, target::Int)
    (pred[target] == 0 && target != source) && return Int[]
    path = Int[target]; curr = target
    while curr != source && curr != 0; curr = pred[curr]; curr == 0 && return Int[]
        pushfirst!(path, curr)
    end
    return path
end
end
