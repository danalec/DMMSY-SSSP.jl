module CSRGraphModule

export CSRGraph, reconstruct_path, random_graph, export_dot

"""
    CSRGraph{T<:Integer, W<:Real}

Compressed Sparse Row representation of a directed graph.
"""
struct CSRGraph{T<:Integer, W<:Real}
    n::T                      # Number of vertices
    m::T                      # Number of edges
    offset::Vector{T}         # offset[i] = start index of adjacency list for vertex i
    adjacency::Vector{T}      # Target vertices (flattened adjacency lists)
    weights::Vector{W}        # Edge weights (parallel to adjacency)
    mean_weight::W            # Pre-computed mean edge weight
end

function CSRGraph{T, W}(n::Int, edges::Vector{Tuple{Int, Int, W}}) where {T<:Integer, W<:Real}
    buckets = [Tuple{Int, W}[] for _ in 1:n]
    for (u, v, w) in edges
        if 1 <= u <= n && 1 <= v <= n && w >= 0
            push!(buckets[u], (v, w))
        end
    end

    m = sum(length(bucket) for bucket in buckets)
    offset = zeros(T, n + 1)
    adjacency = zeros(T, m)
    weights = zeros(W, m)

    idx = 1
    for u in 1:n
        offset[u] = idx
        for (v, w) in buckets[u]
            adjacency[idx] = v
            weights[idx] = w
            idx += 1
        end
    end
    offset[n + 1] = m + 1

    mean_w = m > 0 ? sum(weights) / m : zero(W)
    return CSRGraph{T, W}(T(n), T(m), offset, adjacency, weights, mean_w)
end

CSRGraph(n::Int, edges::Vector{Tuple{Int, Int, Float64}}) = CSRGraph{Int, Float64}(n, edges)
CSRGraph(n::Int, edges::Vector{Tuple{Int, Int, W}}) where W<:Real = CSRGraph{Int, W}(n, edges)

"""
    reconstruct_path(pred::Vector{Int}, source::Int, target::Int)

Reconstruct path from source to target using predecessor array.
"""
function reconstruct_path(pred::Vector{T}, source::T, target::T) where T<:Integer
    (pred[target] == 0 && target != source) && return T[]
    path = T[target]
    current = target
    while current != source && current != 0
        current = pred[current]
        current == 0 && return T[]
        pushfirst!(path, current)
    end
    return path
end

"""
    random_graph(n::Int, m::Int, max_weight::Real)

Generate a random directed graph.
"""
function random_graph(n::Int, m::Int, max_weight::W) where W<:Real
    edges = Tuple{Int, Int, W}[]
    if n <= 0
        return CSRGraph{Int, W}(0, edges)
    end
    for _ in 1:m
        push!(edges, (rand(1:n), rand(1:n), rand() * max_weight))
    end
    return CSRGraph{Int, W}(n, edges)
end

"""
    export_dot(g::CSRGraph, filename::String="graph.dot")

Export the graph to a Graphviz DOT file for visualization.
"""
function export_dot(g::CSRGraph{T, W}, filename::String="graph.dot") where {T, W}
    open(filename, "w") do io
        println(io, "digraph G {")
        println(io, "  node [shape=circle];")
        for u in 1:g.n
            start = g.offset[u]
            stop = g.offset[u+1] - 1
            for i in start:stop
                v = g.adjacency[i]
                w = g.weights[i]
                println(io, "  $u -> $v [label=\"$(round(w, digits=2))\"];")
            end
        end
        println(io, "}")
    end
    println("Graph exported to $filename")
end

end # module
