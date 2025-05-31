module Suurballe

using Graphs
using SimpleWeightedGraphs

struct SuurballeSettings
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

function suurballe(g, runS::SuurballeSettings)
    edges_mat = adjacency_matrix(g)
    hits = ones(size(edges_mat)) .* 0.0000001

    # Adds single edge CFP-s as weights to the graph as weights
    for (break_val, edges) in zip(runS.break_probs, runS.break_edges)
        if length(edges) == 1
            s, e = edges[1]
            hits[s, e] = -log(1 - break_val)
            hits[e, s] = -log(1 - break_val)
        end
    end

    s = runS.start_p
    t = runS.end_p
    ϵ = 0.0001 # make real edges non-zero

    adj_mat = hits .* edges_mat
    l, _ = size(adj_mat)
    g_1 = SimpleWeightedDiGraph(adj_mat)

    # 1. Find the shortest path tree T rooted at node s by running Dijkstra's algorithm
    d = dijkstra_shortest_paths(g_1, s)
    p_1 = vec_from_parents(t, d.parents)
    edges_1 = collect(zip(p_1, p_1[2:end]))

    # 2. Modify the cost of each edge in the graph
    # Type has to be set, otherwise dijsktra does not work since typemax(::Type{Real}) doesn't work
    adj_mat_2::Matrix{Float64} = [adj_mat[i, j] > 0 ? (adj_mat[i, j] - d.dists[j] + d.dists[i]) + ϵ : 0 for i in 1:l, j in 1:l]

    # Create a residual graph G_t formed from G by removing the edges of G on path P1 that are directed into s and then reverse the direction of the zero length edges along path P1.
    i = t
    while d.parents[i] != 0
        j = d.parents[i]
        adj_mat_2[i, j] = ϵ
        adj_mat_2[j, i] = 0
        i = j
    end

    g_t = SimpleWeightedDiGraph(adj_mat_2)

    # 4. Find the shortest path P2 in the residual graph Gt by running Dijkstra's algorithm
    d_2 = dijkstra_shortest_paths(g_t, s)
    p_2 = vec_from_parents(t, d_2.parents)
    edges_2 = collect(zip(p_2, p_2[2:end]))
    # 5. Discard the reversed edges of P2 from both paths. The remaining edges of P1 and P2 form a subgraph with two outgoing edges at s, two incoming edges at t, and one incoming and one outgoing edge at each remaining vertex. Therefore, this subgraph consists of two edge-disjoint paths from s to t and possibly some additional (zero-length) cycles. Return the two disjoint paths from the subgraph.

    adj_mat_3 = [
        (in((i, j), edges_1) && !in((j, i), edges_2)) ||
        (in((i, j), edges_2) && !in((j, i), edges_1)) ? 1 : 0
        for i in 1:l, j in 1:l
    ]
    g_3 = SimpleWeightedDiGraph(adj_mat_3)
    d3 = dijkstra_shortest_paths(g_3, s)
    p_1_final = vec_from_parents(t, d3.parents)
    adj_mat_4 = deepcopy(adj_mat_3)
    for (i, j) in zip(p_1_final, p_1_final[2:end])
        adj_mat_4[i, j] = 0
    end
    g_4 = SimpleWeightedDiGraph(adj_mat_4)
    d_4 = dijkstra_shortest_paths(g_4, s)
    p_2_final = vec_from_parents(t, d_4.parents)

    if length(p_1_final) <= 1
        return (p_2_final, p_2_final)
    end

    if length(p_2_final) <= 1
        return (p_1_final, p_1_final)
    end

    p_1_final, p_2_final
end

end