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
	using Statistics

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
function plot_graph_sets(nodes, edges, fp_edges, solution, title="", color=:black, linestyle=:solid)

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
    p = plot(grid=false, legend=false, aspect_ratio=1.25, axis=false, title=title)

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
    g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/italy_6/italy6.gml", "graphs/italy_6/cfp.xml", "graphs/italy_6/fp.xml")

    (nodes, edges) = read_graph_with_positions("graphs/italy_6/italy6.gml")
    # Example graph data

    ga_better_end = [2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]

    ga_worse_start = [1, 0, 1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
    ga_worse_end = [0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0]
end

# ╔═╡ ab5e803b-9a98-4a67-948c-4287afc50b9c
naive = [1, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 0, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 0, 0, 2, 1, 1, 1, 1, 0, 1, 1, 1, 1, 2, 0, 0, 0, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 0, 1, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 0]

# ╔═╡ 90a5b1b0-33e8-4ecf-8adf-257ae691c0f6
fp_edges[length.(fp_edges).==2]

# ╔═╡ 78940c75-3641-4015-9dc8-ddf26ec83af7
@bind setting sets(fp_edges, ga_worse_start
)

# ╔═╡ 68e3923e-134d-456d-a06d-5d180ac373e7
solution = [getproperty(setting, Symbol("$(i)")) for (i, edges) in enumerate(fp_edges)]

# ╔═╡ 70a858e1-cf99-46bd-8d73-a4f8072bcb96
p = plot_graph_sets(nodes, edges, fp_edges, solution)

# ╔═╡ 5169561b-9f52-4f11-a370-6346d8756a49
sum(fps[solution.!=0])

# ╔═╡ 1958ddf0-bb41-4ea5-a000-e0a101c30dc0
manual = [1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]

# ╔═╡ 98066b1b-5021-4a92-a8c2-e9e53e3bbbd3
italy_5_18_generations = [1, 10, 13, 16, 20, 23, 28, 31, 32, 34]#[1, 5, 10, 13, 14, 34, 39, 41, 46, 48, 51, 52]

# ╔═╡ 66c49583-c6c8-466b-b1d1-ab706dbbc240
italy_5_18_sets = [
[1, 0, 1, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
[0, 0, 1, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]

]


# ╔═╡ e4ca0f51-6374-48f6-9119-3eb6412dca3b
italy_5_18_fitness = [0.12141871101851757,
0.12141898397573535,
0.12143067123385899,
0.12145509327510158,
0.12145707372215196,
0.12362247961010438,
0.12364012817935872,
0.12364017339920218,
0.12443529635431437,
0.12443559934422802]

# ╔═╡ d897679d-b9d2-4512-a055-6a172317a555
@bind version_i Slider(1:length(italy_5_18_generations), show_value=true)

# ╔═╡ 7eeef01c-f9ca-48af-a2ea-9070af5c7d74
plot_graph_sets(nodes, edges, fp_edges, italy_5_18_sets[version_i], "1 -> 20 GEN: $(italy_5_18_generations[version_i]) FIT: $(italy_5_18_fitness[version_i])")

# ╔═╡ 836c363f-ed3f-4405-8952-749001170026
for i in 1:length(italy_5_18_generations)
	savefig(plot_graph_sets(nodes, edges, fp_edges, italy_5_18_sets[i], "5 -> 18 GEN: $(italy_5_18_generations[i]) FIT: $(italy_5_18_fitness[i])"), "italy6_1_20_$(i).pdf")
end

# ╔═╡ de1822f2-ea8f-40b0-a34f-19bd97b66161
function get_line(p1, p2)
	x1, y1 = p1
	x2, y2 = p2
	m = (y2-y1) / (x2-x1)
	n = y1 - m * x1

	m, n
end

# ╔═╡ 6f84650b-18e7-4129-9095-863468f63737
function get_position_compared_to_line(nodes, edges, line)
	if length(edges) > 1
		return 0
	end

	edge = edges[1]
	m, n = line

	s, t = edge

	
	s_side = sign(m * lon(nodes[s]) + n - lat(nodes[s]))
	t_side = sign(m * lon(nodes[t]) + n - lat(nodes[t]))
	if s_side != t_side
		if (s_side == 0 && t_side > 0) || (t_side == 0 && s_side > 0)
			return 1
		elseif (s_side == 0 && t_side < 0) || (t_side == 0 && s_side < 0)
			return 2
		end
		return 0
	elseif s_side > 0
		return 1
	end

	return 2
end

# ╔═╡ 6735e6fb-e5a2-4b53-b290-8a31558e44f2
@bind line_s Slider(1:length(nodes), show_value=true)

# ╔═╡ 5f36e65c-60c6-480d-8d63-4f59a67820c9
@bind line_t Slider(1:length(nodes), show_value=true)

# ╔═╡ a0f8d6fb-9768-41e3-be1d-d6e7d9548838
line = get_line((lon(nodes[line_s]), lat(nodes[line_s])), (lon(nodes[line_t]), lat(nodes[line_t])))

# ╔═╡ 918c6c82-8df8-4a80-8570-7921ad41b4c6
line_sets = [get_position_compared_to_line(nodes, x, line) for x in fp_edges]

# ╔═╡ 57601913-9619-401e-90ec-d3fe4dc18991
begin
	p_line = plot_graph_sets(nodes, edges, fp_edges, line_sets, "$(line_s) -> $(line_t)")
	plot!(p_line, [lon(nodes[line_s]), lon(nodes[line_t])], [lat(nodes[line_s]), lat(nodes[line_t])], color=:green, width=4)
end

# ╔═╡ 95bc7711-86eb-4f6f-b0a6-7b3fc0af130c
runS1 = SafestRoutePair.GaRunSettings(g, fps, fp_edges, [], cfps, cfp_edges, line_s, line_t)

# ╔═╡ a46be0a8-3760-414e-bee0-1a1a519e585c
gaS = GeneticSettings(50, 0.5, 0.9, 0.5, npoint_crossover_naive, mutate_random_multiple, 100, calc_fitness_sets, 0, 0, false)

# ╔═╡ d55a3176-2c69-455c-9a21-8b547f987f5c
chromosomes = [ i < 3 ? line_sets : [length(es) == 1 ? rand([0, 1, 2]) : 0 for es in runS1.fp_edges] for i in 1:gaS.populationSize] #new_fp_edges

# ╔═╡ 0d7ff4d9-9c89-4e60-a543-70e41d3012ee
res = [SafestRoutePair.genetic(runS1, gaS, chromosomes; logging_file="logs_$(i).csv", use_folds=true) for i in 1:50]

# ╔═╡ 045b1b67-e20d-4410-8636-a50f5fd4b3d0
function calc_edges_from_nodes(x)
    collect(zip(x, x[2:end]))
end

# ╔═╡ d101725c-f3aa-4add-94f2-683fbf978bd3
fst((x, _)) = x

# ╔═╡ 153db1a8-0efb-4287-b7fd-d66928b181a7
snd((_, y)) = y

# ╔═╡ b290a757-61ef-436c-a219-92a4d91491f5
function path_intersects(xs, ys)
    for y in ys

        if findfirst(fst.(xs) .== fst(y) .&& snd.(xs) .== snd(y)) !== nothing ||
           findfirst(fst.(xs) .== snd(y) .&& snd.(xs) .== fst(y)) !== nothing
            return true
        end
    end

    false
end

# ╔═╡ fe2fb647-e315-4b2f-83f1-a6a2dac9af82
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

# ╔═╡ 0cefbbd8-1e50-4499-872c-434c1704ca16
res_vals = [calc_availability(x, fps, fp_edges) for x in res]

# ╔═╡ b5d929a7-dfc3-4213-9ef6-e36c6eca921c
mean(res_vals), std(res_vals), maximum(res_vals), minimum(res_vals)

# ╔═╡ f58ef333-1442-4468-b893-8fa267e34eb5
res_naive = SafestRoutePair.safest_route_pair_naive(SafestRoutePair.GraphWithFPandCFP(g, fps, fp_edges, cfps, cfp_edges), line_s, line_t)

# ╔═╡ c87e2050-c14d-48f4-a9c6-b9dbf296df50
res_val_naive = calc_availability(fst(res_naive), fps, fp_edges)

# ╔═╡ a0c7be8d-ae05-4a27-a890-7664e7c43c63
res_suurballe = SafestRoutePair.safest_route_pair_suurballe(SafestRoutePair.GraphWithFPandCFP(g, fps, fp_edges, cfps, cfp_edges), line_s, line_t)

# ╔═╡ d95d4923-1f37-4e49-9e6e-09b6ee6dc5e8
res_val_suurballe = calc_availability(fst(res_suurballe), fps, fp_edges)

# ╔═╡ 1072dc58-b726-40d6-9c0b-dfd5b67e8748
res_aco = [ SafestRoutePair.safest_route_pair_aco(SafestRoutePair.GraphWithFPandCFP(g, fps, fp_edges, cfps, cfp_edges), line_s, line_t, SafestRoutePair.ACOSettings(1.5, 0.5, 25, 0.1, 0.3, 50, 1)) for i in 1:50]

# ╔═╡ 6554d132-4af1-4cbf-9780-3899477379de
res_vals_aco = [calc_availability(fst(x), fps, fp_edges) for x in res_aco]

# ╔═╡ 8d8d1799-465b-4b15-bb9d-d9f6e523369d
mean(res_vals_aco), std(res_vals_aco), maximum(res_vals_aco), minimum(res_vals_aco)

# ╔═╡ Cell order:
# ╠═606e0d46-12c5-11f0-31be-7b0572b32ee3
# ╠═ff68e30f-66f1-4cc6-a887-312f24390492
# ╠═1c31ef47-7995-4a1b-83d6-93331d819062
# ╠═7d17dc41-2e03-48d6-89ae-5f367b9377c8
# ╠═ab5e803b-9a98-4a67-948c-4287afc50b9c
# ╠═90a5b1b0-33e8-4ecf-8adf-257ae691c0f6
# ╠═78940c75-3641-4015-9dc8-ddf26ec83af7
# ╠═68e3923e-134d-456d-a06d-5d180ac373e7
# ╠═70a858e1-cf99-46bd-8d73-a4f8072bcb96
# ╠═5169561b-9f52-4f11-a370-6346d8756a49
# ╠═1958ddf0-bb41-4ea5-a000-e0a101c30dc0
# ╠═98066b1b-5021-4a92-a8c2-e9e53e3bbbd3
# ╠═66c49583-c6c8-466b-b1d1-ab706dbbc240
# ╠═e4ca0f51-6374-48f6-9119-3eb6412dca3b
# ╠═d897679d-b9d2-4512-a055-6a172317a555
# ╠═7eeef01c-f9ca-48af-a2ea-9070af5c7d74
# ╠═836c363f-ed3f-4405-8952-749001170026
# ╠═de1822f2-ea8f-40b0-a34f-19bd97b66161
# ╠═6f84650b-18e7-4129-9095-863468f63737
# ╠═6735e6fb-e5a2-4b53-b290-8a31558e44f2
# ╠═5f36e65c-60c6-480d-8d63-4f59a67820c9
# ╠═a0f8d6fb-9768-41e3-be1d-d6e7d9548838
# ╠═918c6c82-8df8-4a80-8570-7921ad41b4c6
# ╠═57601913-9619-401e-90ec-d3fe4dc18991
# ╠═95bc7711-86eb-4f6f-b0a6-7b3fc0af130c
# ╠═a46be0a8-3760-414e-bee0-1a1a519e585c
# ╠═d55a3176-2c69-455c-9a21-8b547f987f5c
# ╠═0d7ff4d9-9c89-4e60-a543-70e41d3012ee
# ╟─045b1b67-e20d-4410-8636-a50f5fd4b3d0
# ╟─d101725c-f3aa-4add-94f2-683fbf978bd3
# ╟─153db1a8-0efb-4287-b7fd-d66928b181a7
# ╟─b290a757-61ef-436c-a219-92a4d91491f5
# ╟─fe2fb647-e315-4b2f-83f1-a6a2dac9af82
# ╠═0cefbbd8-1e50-4499-872c-434c1704ca16
# ╠═b5d929a7-dfc3-4213-9ef6-e36c6eca921c
# ╠═f58ef333-1442-4468-b893-8fa267e34eb5
# ╠═c87e2050-c14d-48f4-a9c6-b9dbf296df50
# ╠═a0c7be8d-ae05-4a27-a890-7664e7c43c63
# ╠═d95d4923-1f37-4e49-9e6e-09b6ee6dc5e8
# ╠═1072dc58-b726-40d6-9c0b-dfd5b67e8748
# ╠═6554d132-4af1-4cbf-9780-3899477379de
# ╠═8d8d1799-465b-4b15-bb9d-d9f6e523369d
