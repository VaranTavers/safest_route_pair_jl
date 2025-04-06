using Pkg

Pkg.activate(".")

using Plots
using Colors

using SafestRoutePair

import SafestRoutePair: read_graph_with_positions


#=
nodes, edges = read_graph_with_positions("graphs/italy7/italy7.gml")

lines = graph_to_tikz_net(nodes, edges)
for l in lines
    println(l)
end
=#


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

function plot_graph_set_minus(nodes, edges, fp_edges, solution, color=:black, linestyle=:solid)

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

g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/italy7/italy7.gml", "graphs/italy7/cfp.xml", "graphs/italy7/fp.xml")
(nodes, edges) = read_graph_with_positions("graphs/italy7/italy7.gml")
# Example graph data

ga_better_end = [2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]

ga_worse_start = [1, 0, 1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
ga_worse_end = [0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0]
p = plot_graph_sets(nodes, edges, fp_edges, ga_worse_end)
# Display plot
savefig(p, "test_4.pdf")