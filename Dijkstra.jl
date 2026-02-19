module DijkstraModule

using DataStructures: BinaryMinHeap
using ..CSRGraphModule

export dijkstra_ref

"""
    dijkstra_ref(g::CSRGraph, source::Integer)

Standard Dijkstra implementation for reference.
"""
function dijkstra_ref(g::CSRGraph{T, W}, source::T) where {T<:Integer, W<:Real}
    n = g.n
    d = fill(typemax(W), n)
    pred = zeros(T, n)
    visited = falses(n)

    d[source] = zero(W)
    heap = BinaryMinHeap{Tuple{W, T}}()
    push!(heap, (zero(W), source))
    
    adj = g.adjacency
    weights = g.weights
    offsets = g.offset

    while !isempty(heap)
        (dist_u, u) = pop!(heap)
        visited[u] && continue
        visited[u] = true
        
        start = offsets[u]
        stop = offsets[u+1] - 1

        @inbounds for i in start:stop
            v = adj[i]
            w = weights[i]
            if !visited[v] && d[u] + w < d[v]
                d[v] = d[u] + w
                pred[v] = u
                push!(heap, (d[v], v))
            end
        end
    end

    return d, pred
end

end # module
