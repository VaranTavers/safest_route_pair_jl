module TIKZPlotter

using ..GraphReader

import ..GraphReader: GeoPoint, NodeWithGeoPoint, EdgeWithGeoPoints

function vertex_str(n::NodeWithGeoPoint; label=true, id=false)

    @assert !(label && id) "Label and id should not be used together"

    r_label = label ? n.label : (id ? n.id : "")
    r_label = replace(r_label, '_' => "\\_")
    "\\Vertex[y=$(n.point.lat), x=$(n.point.lon), label=$(r_label)]{$(n.id)}"
end

function node_str(p::GeoPoint, label; draw=false)
    "\\Vertex[y=$(p.lat), x=$(p.lon), shape = coordinate]{$(label)}"
    #"\\node[x=$(p.lat), y=$(p.lon), $(draw ? "" : "draw = none")]($(label)){};"
end

function get_point(p::GeoPoint)
    "{$(p.lon), $(p.lat)}"
end

function get_points(ps::Vector{GeoPoint})
    join(get_point.(ps), ',')
end

function edge_str(e::EdgeWithGeoPoints; id=false, bend=0, style="solid", color="black")::String
    if length(e.points) == 0
        return "\\Edge[$(id ? "label = $(e.id), " : "")bend = $(bend), style={$(style)}, color=$(color)]($(e.source))($(e.target))"
    end

    "\\Edge[path={$(e.source),$(get_points(e.points)),$(e.target)}, $(id ? "label = $(e.id), " : "")bend = $(bend), style={$(style)}, color=$(color)]($(e.source))($(e.target))"
end

function graph_to_tikz_net(nodes, edges)
    res = []

    for n in nodes
        push!(res, vertex_str(n; label=true, id=false))
    end

    for e in edges
        push!(res, edge_str(e)) # ; style="dashed", color="red"
    end

    res
end
end