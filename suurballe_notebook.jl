### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ e5451ea6-3d72-11f0-2712-cf3c27927db0
begin
	using Pkg

    Pkg.activate(".")
	using SafestRoutePair
	using Graphs
	using SimpleWeightedGraphs
end

# ╔═╡ 5756790e-040e-44d1-8a44-898eaa2de36d
g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/italy_6/italy6.gml", "graphs/italy_6/cfp.xml", "graphs/italy_6/fp.xml")

# ╔═╡ 2020f7f4-6540-40a2-945c-4d0ad365f412
begin
	edges_mat = adjacency_matrix(g)
    hits = ones(size(edges_mat)) .* 0.0000001

    # Adds single edge CFP-s as weights to the graph as weights
    for (break_val, edges) in zip(cfps, cfp_edges)
        if length(edges) == 1
            s, e = edges[1]
            hits[s, e] = -log(1-break_val)
            hits[e, s] = -log(1-break_val)
        end
    end
end

# ╔═╡ 6d24686d-f547-41b0-b37b-3d6e4a95bd44
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

# ╔═╡ e4d6b836-916a-4ae6-b76d-ec7e47c79f08
begin
	s = 7
	t = 24
	ϵ = 0.0001 # make real edges non-zero
end

# ╔═╡ 00a5f9f1-9871-4af7-b558-cb8bd257cb7a
#= example on wiki adj_mat = [
	0 1 2 0 0 0;
	1 0 0 1 2 0;
	2 0 0 2 0 0;
	0 1 2 0 0 1;
	0 2 0 0 0 2;
	0 0 0 1 2 0;
] =#
adj_mat = hits .* edges_mat

# ╔═╡ 9f65b266-a650-4215-8d6f-1b4a12fa34e7
l, _ = size(adj_mat)

# ╔═╡ 4546dc8b-8fe2-4224-9f57-ea20c1080cdc
g_1 = SimpleWeightedDiGraph(adj_mat)

# ╔═╡ b68dc771-61ac-4061-9605-ade827e19c81
# 1. Find the shortest path tree T rooted at node s by running Dijkstra's algorithm
d = dijkstra_shortest_paths(g_1, s)

# ╔═╡ 3860309b-68e6-4327-bb7d-91cf6a5c7e92
p_1 = vec_from_parents(t, d.parents)

# ╔═╡ db40dbc5-f21c-47b2-994e-75378ebbcc65
edges_1 = collect(zip(p_1, p_1[2:end]))

# ╔═╡ 4b4c5fca-cee6-4848-ac6d-6b3d1c685165
begin
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

	adj_mat_2
end

# ╔═╡ 71a07c36-64e5-4b21-8503-cfe99b0d4fbf
g_t = SimpleWeightedDiGraph(adj_mat_2)

# ╔═╡ a601aba5-6fb7-4e29-af1c-39ecf154a4e1
# 4. Find the shortest path P2 in the residual graph Gt by running Dijkstra's algorithm
d_2 = dijkstra_shortest_paths(g_t, s)

# ╔═╡ a594a718-1b68-4ef5-a268-79a30173717a
p_2 = vec_from_parents(t, d_2.parents)

# ╔═╡ 8d9dbc75-699f-43b1-998e-a07cdfc3dd19
edges_2 = collect(zip(p_2, p_2[2:end]))

# ╔═╡ 2935982b-08b1-40c1-b587-12745a23281a
# 5. Discard the reversed edges of P2 from both paths. The remaining edges of P1 and P2 form a subgraph with two outgoing edges at s, two incoming edges at t, and one incoming and one outgoing edge at each remaining vertex. Therefore, this subgraph consists of two edge-disjoint paths from s to t and possibly some additional (zero-length) cycles. Return the two disjoint paths from the subgraph.

# ╔═╡ 21a5d654-2100-4d78-ac0b-5ae428f5411e
adj_mat_3 = [
	(in((i, j), edges_1) && !in((j, i), edges_2)) ||
	(in((i, j), edges_2) && !in((j, i), edges_1)) ? 1 : 0
	for i in 1:l, j in 1:l
]

# ╔═╡ 93ea9cb5-d632-40d5-870c-c603ab116f37
g_3 = SimpleWeightedDiGraph(adj_mat_3)

# ╔═╡ cf94192f-de4e-4213-8acf-91489cc58a5f
d3 = dijkstra_shortest_paths(g_3, s)

# ╔═╡ c9fa8c7c-1abc-481c-9cde-44c6830504f5
p_1_final = vec_from_parents(t, d3.parents)

# ╔═╡ 5de785a5-8310-4b4e-8fb7-3ad3a18828a1
begin
	adj_mat_4 = deepcopy(adj_mat_3)
	for (i, j) in zip(p_1_final, p_1_final[2:end])
		adj_mat_4[i, j] = 0
	end
	adj_mat_4
end

# ╔═╡ 2ca6c05d-e7ae-435b-ae2d-8c678fac1474
g_4 = SimpleWeightedDiGraph(adj_mat_4)

# ╔═╡ 8da8df87-ed79-4993-9cce-2ca8d7f3a62d
d_4 = dijkstra_shortest_paths(g_4, s)

# ╔═╡ 2b92f68f-d586-4884-be40-51cb86c00ad4
p_2_final = vec_from_parents(t, d_4.parents)

# ╔═╡ Cell order:
# ╠═e5451ea6-3d72-11f0-2712-cf3c27927db0
# ╠═5756790e-040e-44d1-8a44-898eaa2de36d
# ╠═2020f7f4-6540-40a2-945c-4d0ad365f412
# ╠═6d24686d-f547-41b0-b37b-3d6e4a95bd44
# ╠═e4d6b836-916a-4ae6-b76d-ec7e47c79f08
# ╠═00a5f9f1-9871-4af7-b558-cb8bd257cb7a
# ╠═9f65b266-a650-4215-8d6f-1b4a12fa34e7
# ╠═4546dc8b-8fe2-4224-9f57-ea20c1080cdc
# ╠═b68dc771-61ac-4061-9605-ade827e19c81
# ╠═3860309b-68e6-4327-bb7d-91cf6a5c7e92
# ╠═db40dbc5-f21c-47b2-994e-75378ebbcc65
# ╠═4b4c5fca-cee6-4848-ac6d-6b3d1c685165
# ╠═71a07c36-64e5-4b21-8503-cfe99b0d4fbf
# ╠═a601aba5-6fb7-4e29-af1c-39ecf154a4e1
# ╠═a594a718-1b68-4ef5-a268-79a30173717a
# ╠═8d9dbc75-699f-43b1-998e-a07cdfc3dd19
# ╠═2935982b-08b1-40c1-b587-12745a23281a
# ╠═21a5d654-2100-4d78-ac0b-5ae428f5411e
# ╠═93ea9cb5-d632-40d5-870c-c603ab116f37
# ╠═cf94192f-de4e-4213-8acf-91489cc58a5f
# ╠═c9fa8c7c-1abc-481c-9cde-44c6830504f5
# ╠═5de785a5-8310-4b4e-8fb7-3ad3a18828a1
# ╠═2ca6c05d-e7ae-435b-ae2d-8c678fac1474
# ╠═8da8df87-ed79-4993-9cce-2ca8d7f3a62d
# ╠═2b92f68f-d586-4884-be40-51cb86c00ad4
