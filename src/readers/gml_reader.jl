module GMLReader

using Graphs
using SimpleWeightedGraphs
using GraphIO
using ParserCombinator

import Graphs: loadgraph, loadgraphs, savegraph

export CustomGMLFormat

struct CustomGMLFormat <: Graphs.AbstractGraphFormat
    lists::Any
    CustomGMLFormat() = new([:graph, :node, :edge, :point])
end

function _gml_read_one_graph(gs, dir)
    nodes = [x[:id] for x in gs[:node]]
    if dir
        g = DiGraph(length(nodes))
    else
        g = Graph(length(nodes))
    end
    mapping = Dict{Int,Int}()
    # If the original nodes are not in the right order it could mess up the pairings of routes and their CFP-s 
    sort!(nodes)
    #@show nodes
    for (i, n) in enumerate(nodes)
        mapping[n] = i
    end
    #@show mapping
    sds = [(Int(x[:source]), Int(x[:target])) for x in gs[:edge]]
    for (s, d) in (sds)
        #@show s, d, mapping[s], mapping[d]
        add_edge!(g, mapping[s], mapping[d])
    end

    return g
end

function loadgml(io::IO, gname::String, lists)
    p = Parsers.GML.parse_dict(read(io, String), lists=lists)
    for gs in p[:graph]
        dir = Bool(get(gs, :directed, 0))
        graphname = get(gs, :label, dir ? "digraph" : "graph")

        (gname == graphname) && return _gml_read_one_graph(gs, dir)
    end
    return error("Graph $gname not found")
end

loadgraph(io::IO, gname::String, format::CustomGMLFormat) = loadgml(io, gname, format.lists)


# Custom implementation for graph plotting

function parseGMLDict(f; depth=0, log=false)
    ret = Dict()
    row = readline(f)

    while !(occursin("]", row))
        if isempty(row) && depth == 0
            break
        end
        if !startswith(row, repeat('\t', depth)) && log
            println("Bad indentation on line: $(row)")
        end
        if !isempty(row) && !isempty(strip(row))
            parts = split(strip(row), ' ', limit=2)
            propName = parts[1]
            if propName in keys(ret) && !(ret[propName] isa AbstractVector)
                ret[propName] = [ret[propName]]
            end
            if occursin("[", row)
                propVal = parseGMLDict(f; depth=depth + 1, log=log)
            else
                if occursin("\"", row)
                    propVal = replace(parts[2], '"' => "")
                elseif occursin(".", row)
                    propVal = parse(Float64, parts[2])
                else
                    propVal = parse(Int64, parts[2])
                end
            end
            if propName in keys(ret)
                push!(ret[propName], propVal)
            else
                ret[propName] = propVal
            end
        end

        row = readline(f)
    end
    ret
end


struct GeoPoint
    lon::Float64
    lat::Float64
end

struct NodeWithGeoPoint
    point::GeoPoint
    label::String
    id::Int64
end

struct EdgeWithGeoPoints
    source::Int64
    target::Int64
    id::Int64
    points::Vector{GeoPoint}
end

function read_points(points_dict)
    ret = []

    for point in points_dict["point"]
        push!(ret, GeoPoint(point["Longitude"][1], point["Latitude"][1]))
    end

    ret
end

function read_longitude(node_dict)
    possible_names = ["Longitude", "longitude", "lon"]
    for name in possible_names
        if name in keys(node_dict)
            return node_dict[name]
        end
    end

    @error("Latitude information missing from node $(node_dict)")
end

function read_latitude(node_dict)
    possible_names = ["Latitude", "latitude", "lat"]
    for name in possible_names
        if name in keys(node_dict)
            return node_dict[name]
        end
    end

    @error("Latitude information missing from node $(node_dict)")
end

function create_node(node_dict)

    NodeWithGeoPoint(
        GeoPoint(
            read_longitude(node_dict),
            read_latitude(node_dict)
        ),
        node_dict["label"],
        node_dict["id"]
    )
end

function create_edge(edge_dict, i=1)
    points = []
    id = i
    if "points" in keys(edge_dict)
        points = read_points(edge_dict["points"])
    end
    if "id" in keys(edge_dict)
        id = edge_dict["id"]
    end

    EdgeWithGeoPoints(
        edge_dict["source"],
        edge_dict["target"],
        id,
        points
    )
end

function read_graph_with_positions(file; log=false)
    nodes = []
    edges = []
    open(file) do f
        dict = parseGMLDict(f; log=log)["graph"]
        for node in dict["node"]
            push!(nodes, create_node(node))
        end
        for (i, edge) in enumerate(dict["edge"])
            push!(edges, create_edge(edge, i))
        end
    end

    sort!(nodes, lt=(x, y) -> x.id < y.id)

    nodes, edges
end


lat(node::NodeWithGeoPoint) = node.point.lat
lon(node::NodeWithGeoPoint) = node.point.lon

lat(point::GeoPoint) = point.lat
lon(point::GeoPoint) = point.lon


end