### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 82c977da-583a-11f0-3fed-fbc987915aca
begin 
	# Loads the local package and the required dependencies.
	using Pkg
	Pkg.activate(".")

	using SafestRoutePair
	using Plots
	using Colors
	using Statistics
	using Random
	using PlutoUI
	using DataFrames
end

# ╔═╡ 5b839d20-95c5-44a9-89c0-3352c2bae7e4
md"""
# Safest route pair

Finds a pair of routes between the source and the target that has maximal uptime
"""

# ╔═╡ 8142fd0c-18c9-4d10-99a4-4668859d2646
begin
	@assert isdir("./graphs") "Graphs folder does not exist"

	fst((x,_),) = x
	snd((_,y),) = y
	
	function calc_edges_from_nodes(x)
    	collect(zip(x, x[2:end]))
	end

	function path_intersects(xs, ys)
    	for y in ys
        	if findfirst(fst.(xs) .== fst(y) .&& snd.(xs) .== snd(y)) !== nothing ||
           	findfirst(fst.(xs) .== snd(y) .&& snd.(xs) .== fst(y)) !== nothing
            	return true
        	end
    	end

    	false
	end
	
	function calc_availability((path1, path2), fps, fp_edges)
    # Transform solution from list of nodes into list of edges (node pairs)
    path_a = calc_edges_from_nodes(path1)
    path_b = calc_edges_from_nodes(path2)

    only_path_a = 0
    only_path_b = 0
    both_paths = 0

    for (prob, vec) in zip(fps, fp_edges)
        intersects_path_a = path_intersects(path_a, vec)
        intersects_path_b = path_intersects(path_b, vec)
        if intersects_path_a && intersects_path_b
            both_paths += prob
        elseif intersects_path_a
            only_path_a += prob
        elseif intersects_path_b
            only_path_b += prob
        end
    end


    res = both_paths + only_path_a * only_path_b


    -log(res)
end
	
	# Drawing functions
	function plot_graph_routes(nodes, edges, route_1, route_2, use_labels, title, linestyle=:solid)
		n = length(nodes)
		node_pairs_1 = collect(zip(route_1, route_1[2:end]))
		node_pairs_2 = collect(zip(route_2, route_2[2:end]))
	
    	# Initialize plot
    	p = plot(grid=false, legend=false, aspect_ratio=1.25, axis=false, title=title)

    	# Plot edges
    	for edge in edges
        	lats = [lat(p) for p in edge.points]
        	lons = [lon(p) for p in edge.points]

        	source = edge.source + 1
        	target = edge.target + 1

			pushfirst!(lats, lat(nodes[source]))
			push!(lats, lat(nodes[target]))
			pushfirst!(lons, lon(nodes[source]))
			push!(lons, lon(nodes[target]))

			in_p1 = (source, target) in node_pairs_1  || (target, source) in node_pairs_1
			in_p2 = (target, source) in node_pairs_2 || (source, target) in node_pairs_2
			
        	if in_p1 && in_p2
            	plot!(p, lons, lats, label="Both paths", linewidth=2, axis=false, color=:purple)
        	elseif in_p1
            	plot!(p, lons, lats, label="Path 1", linewidth=2, axis=false, color=:blue)
			elseif in_p2
            	plot!(p, lons, lats, label="Path 2", linewidth=2, axis=false, color=:red)
			else
            	plot!(p, lons, lats, label="None", linewidth=2, axis=false, color=:gray67)
        	end

    	end

    # Plot nodes
    	scatter!(p, lon.(nodes), lat.(nodes), 
				 markersize=8, 
				 markercolor=[i == route_1[1] || i == route_1[end] ? :red : :lightblue for i in 1:n], 
				 series_annotations=[use_labels ? Plots.text(nodes[i].label, pointsize=6) : Plots.text(nodes[i].id, pointsize=6) for i in 1:n], 
				 axis=false
				)

    	p
	end
	
	graphs = readdir("./graphs")
end

# ╔═╡ 726509e8-8f53-4ce5-b791-ed9eb1b0f955
md"Please select a graph: $(@bind graph Select(graphs))"

# ╔═╡ 0ca334e1-fa07-4827-b5e6-1c6914f2b4fd
md"Use node labels $(@bind node_labels CheckBox(default=false))"

# ╔═╡ 1c6e2e42-dd5b-44f0-a035-dcbfbfc29d39
begin
	files = readdir("./graphs/$(graph)")
	gml_file = files[findfirst(x -> x[end-3:end] == ".gml", files)]
	(nodes, edges) = read_graph_with_positions("graphs/$(graph)/$(gml_file)")
	g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/$(graph)/$(gml_file)", "graphs/$(graph)/cfp.xml", "graphs/$(graph)/fp.xml")
	gwfp = GraphWithFPandCFP(g, fps, fp_edges, cfps, cfp_edges)
end;

# ╔═╡ 00a03837-9ba9-46d9-9465-699b06eae046
begin
	node_id_plus_label = sort([(nodes[i].id + 1) => "$(nodes[i].id):$(nodes[i].label)" for i in 1:length(nodes)])
	md"From: $(@bind node_s Select(node_id_plus_label)) To: $(@bind node_t Select(node_id_plus_label))" 
end

# ╔═╡ c1c94bdb-4a82-4f09-8ba0-5a00eca8a707
md"α: $(@bind α Slider(0.5:0.1:2, default=1.5, show_value=true)) β: $(@bind β Slider(0.5:0.1:2, show_value=true)) ants: $(@bind ants Slider(10:5:100, show_value=true)) gens: $(@bind gens Slider(10:5:100, show_value=true))"

# ╔═╡ 756ead37-c75b-410e-a72f-02b2c069a1f1
begin
	@show "Naive runtime:"
	naive, _ = @time SafestRoutePair.safest_route_pair_naive(gwfp, node_s, node_t)
	@show "Suurballe runtime:"
	suurballe, _ = @time SafestRoutePair.safest_route_pair_suurballe(gwfp, node_s, node_t)
	@show "ACO runtime:"
	aco, _ = @time SafestRoutePair.safest_route_pair_aco(gwfp, node_s,  node_t, ACOSettings(α, β, ants, 0.3, 0.1, gens, 1))
	routes = Dict("naive" => naive, "suurballe" => suurballe, "ACO" => aco)
end

# ╔═╡ 879e2386-2a9e-45a5-bc2f-dfc83d03926f
begin
	naive_nl = calc_availability(naive, fps, fp_edges)
	suurballe_nl = calc_availability(suurballe, fps, fp_edges)
	aco_nl = calc_availability(aco, fps, fp_edges)
	df = DataFrame("" => ["-log(1-A(P_1, P,2))", "A(P1, P2)"], "naive" => [naive_nl, 1 - exp(-naive_nl)], "suurballe" => [suurballe_nl, 1 - exp(-suurballe_nl)], "ACO" => [aco_nl, 1 - exp(-aco_nl)])
end

# ╔═╡ e1112236-428f-4b71-829d-1b2d4420cd47
md"Paths from algorithm: $(@bind alg Select(collect(keys(routes))))"

# ╔═╡ 84ae099d-b2fc-4fe6-ab8a-7db1ba52c8f3
plot_graph_routes(nodes, edges, fst(routes[alg]), snd(routes[alg]), node_labels, graph)

# ╔═╡ Cell order:
# ╠═5b839d20-95c5-44a9-89c0-3352c2bae7e4
# ╟─82c977da-583a-11f0-3fed-fbc987915aca
# ╟─8142fd0c-18c9-4d10-99a4-4668859d2646
# ╠═726509e8-8f53-4ce5-b791-ed9eb1b0f955
# ╠═0ca334e1-fa07-4827-b5e6-1c6914f2b4fd
# ╟─00a03837-9ba9-46d9-9465-699b06eae046
# ╟─1c6e2e42-dd5b-44f0-a035-dcbfbfc29d39
# ╟─c1c94bdb-4a82-4f09-8ba0-5a00eca8a707
# ╟─756ead37-c75b-410e-a72f-02b2c069a1f1
# ╟─879e2386-2a9e-45a5-bc2f-dfc83d03926f
# ╟─e1112236-428f-4b71-829d-1b2d4420cd47
# ╠═84ae099d-b2fc-4fe6-ab8a-7db1ba52c8f3
