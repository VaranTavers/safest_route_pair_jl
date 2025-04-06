### A Pluto.jl notebook ###
# v0.20.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 606e0d46-12c5-11f0-31be-7b0572b32ee3
begin
	using Pkg

	Pkg.activate(".")

	using Plots
	using SafestRoutePair
	using Plots
	using Colors
	using PlutoUI
	using HypertextLiteral

	import SafestRoutePair: read_graph_with_positions
	import PlutoUI: combine

end

# ╔═╡ ff68e30f-66f1-4cc6-a887-312f24390492
function sets(fp_edges::Vector, defaults::Vector)
	fp_edges = [("$(i)", "$(x)", def) for ((i, x), def) in zip(enumerate(fp_edges), defaults)]
	return combine() do Child
		@htl("""
		<h6>Wind speeds</h6>
		<table>
		$([
			@htl("<tr><td>$(name)</td><td> $(Child(i, Select([0,1,2], default=def)))</td></tr>")
			for (i, name, def) in fp_edges
		])
		</table>
		""")
	end
end

# ╔═╡ 1c31ef47-7995-4a1b-83d6-93331d819062
function plot_graph_sets(nodes, edges, fp_edges, solution, color=:black, linestyle=:solid)

    l_n = length(nodes)
    is_in_s1 = zeros(l_n, l_n)
    is_in_s2 = zeros(l_n, l_n)

    for (sol, fp_edge_list) in zip(solution, fp_edges)
        for (a, b) in fp_edge_list
            if sol == 1
                is_in_s1[a, b] = 1
                is_in_s1[b, a] = 1
            elseif sol == 2

                is_in_s2[a, b] = 1
                is_in_s2[b, a] = 1
            end
        end
    end
    # Initialize plot
    # title="Italy (interroute_v2) $(route[1]) -> $(route[end])"
    p = plot(grid=false, legend=false, aspect_ratio=1.25, axis=false)

    # Plot edges
    for edge in edges
        lats = [lat(p) for p in edge.points]
        lons = [lon(p) for p in edge.points]

        source = edge.source + 1
        target = edge.target + 1

        if is_in_s1[source, target] == 1 && is_in_s2[source, target] == 1
            plot!(p, lons, lats, label="S1 & S2", linewidth=2, axis=false, color=:purple)
        elseif is_in_s1[source, target] == 1
            plot!(p, lons, lats, label="S1", linewidth=2, axis=false, color=:red)
        elseif is_in_s2[source, target] == 1
            plot!(p, lons, lats, label="S2", linewidth=2, axis=false, color=:blue)
        else
            plot!(p, lons, lats, label="None", linewidth=2, axis=false, color=:gray67)
        end

    end

    # Plot nodes
    scatter!(p, lon.(nodes), lat.(nodes), markersize=8, series_annotations=[Plots.text(i - 1, pointsize=6) for i in 1:length(nodes)], axis=false)

    p
end


# ╔═╡ 7d17dc41-2e03-48d6-89ae-5f367b9377c8
begin
g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/italy7/italy7.gml", "graphs/italy7/cfp.xml", "graphs/italy7/fp.xml")

(nodes, edges) = read_graph_with_positions("graphs/italy7/italy7.gml")
# Example graph data

ga_better_end = [2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]

ga_worse_start = [1, 0, 1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
ga_worse_end = [0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0]
end

# ╔═╡ 90a5b1b0-33e8-4ecf-8adf-257ae691c0f6
fp_edges[length.(fp_edges) .== 2]

# ╔═╡ 78940c75-3641-4015-9dc8-ddf26ec83af7
@bind setting sets(fp_edges, [0 for _ in fp_edges])

# ╔═╡ 68e3923e-134d-456d-a06d-5d180ac373e7
solution = [getproperty(setting, Symbol("$(i)")) for (i, edges) in enumerate(fp_edges)]

# ╔═╡ 70a858e1-cf99-46bd-8d73-a4f8072bcb96
p = plot_graph_sets(nodes, edges, fp_edges, solution)

# ╔═╡ 5169561b-9f52-4f11-a370-6346d8756a49
sum(fps[solution .!= 0])

# ╔═╡ 1958ddf0-bb41-4ea5-a000-e0a101c30dc0
manual = [1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]

# ╔═╡ Cell order:
# ╠═606e0d46-12c5-11f0-31be-7b0572b32ee3
# ╠═ff68e30f-66f1-4cc6-a887-312f24390492
# ╠═1c31ef47-7995-4a1b-83d6-93331d819062
# ╠═7d17dc41-2e03-48d6-89ae-5f367b9377c8
# ╠═90a5b1b0-33e8-4ecf-8adf-257ae691c0f6
# ╠═78940c75-3641-4015-9dc8-ddf26ec83af7
# ╠═68e3923e-134d-456d-a06d-5d180ac373e7
# ╠═70a858e1-cf99-46bd-8d73-a4f8072bcb96
# ╠═5169561b-9f52-4f11-a370-6346d8756a49
# ╠═1958ddf0-bb41-4ea5-a000-e0a101c30dc0
