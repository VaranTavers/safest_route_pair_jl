module NaiveComplex

using Graphs
using SimpleWeightedGraphs
# using Folds
# using ConcurrentCollections

using ..GraphUtils
import ..GraphUtils: weighted_graph_from_mat

struct NaiveGreedySettings
    start_p::Integer
    end_p::Integer
    break_probs::Vector{Real}
    break_edges::Vector{Vector{Tuple{Integer,Integer}}}
end


function vec_from_parents(end_p, parents)
    res = zeros(Int64, length(parents))
    res[1] = end_p
    i = 1
    while (parents[res[i]] != 0)
        res[i+1] = parents[res[i]]
        i += 1
    end

    reverse(res[res.!=0])
end

fst((x, _)) = x
snd((_, y)) = y

function are_edges_equal(x, y)
    (fst(x) == fst(y) && snd(x) == snd(y)) || (fst(x) == snd(y) && snd(x) == fst(y))
end

function is_edge_equal_to_any(x, ys)
    for y in ys
        if are_edges_equal(x, y)
            return true
        end
    end

    false
end

function naive_complex(edges_mat, runS::NaiveGreedySettings)
    hits = ones(size(edges_mat)) .* 0.0000001

    # Adds single edge CFP-s as weights to the graph as weights
    for (break_val, edges) in zip(runS.break_probs, runS.break_edges)
        if length(edges) == 1
            s, e = edges[1]
            hits[s, e] = break_val
            hits[e, s] = break_val
        end
    end

    g = weighted_graph_from_mat(hits .* edges_mat)

    D = dijkstra_shortest_paths(g, runS.start_p)

    #@show 1, Matrix(g.weights)

    #@show hits
    #@show D

    path1, dist1 = vec_from_parents(runS.end_p, D.parents), D.dists[runS.end_p]

    edges_path1 = collect(zip(path1, path1[2:end]))

    for (s, e) in edges_path1
        add_edge!(g, s, e, 1 + g.weights[s, e])

    end


    #@show 2, Matrix(g.weights)

    D2 = dijkstra_shortest_paths(g, runS.start_p)

    path2, dist2 = vec_from_parents(runS.end_p, D2.parents), D2.dists[runS.end_p]


    (path1, path2)
end
end