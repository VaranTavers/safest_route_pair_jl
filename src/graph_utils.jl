module GraphUtils

using Graphs
using SimpleWeightedGraphs

function weighted_graph_from_mat(mat)
    n, _ = size(mat)
    g = SimpleWeightedGraph(n)

    for i = 1:n
        for j = 1:n
            add_edge!(g, i, j, mat[i, j])
        end
    end

    g
end

function calc_edges_from_nodes(x)
    collect(zip(x, x[2:end]))
end


function create_indep_graph(g, cfps::Vector{Real}, cfp_edges::Vector{Vector{Tuple{Integer,Integer}}})
    edges_mat = adjacency_matrix(g)
    hits = ones(size(edges_mat)) .* 0.0000001

    # Adds single edge CFP-s as weights to the graph as weights
    for (break_val, edges) in zip(cfps, cfp_edges)
        if length(edges) == 1
            s, e = edges[1]
            hits[s, e] = break_val
            hits[e, s] = break_val
        end
    end

    weighted_graph_from_mat(hits .* edges_mat)
end

end