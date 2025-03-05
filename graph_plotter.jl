using Plots



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
        if !isempty(row)
            parts = split(strip(row), ' ')
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


struct Point
    lon::Float64
    lat::Float64
end

struct Node
    point::Point
    label::String
    id::Int64
end

struct Edge
    source::Int64
    target::Int64
    id::Int64
    points::Vector{Point}
end

function read_points(points_dict)
    ret = []

    for point in points_dict["point"]
        push!(ret, Point(point["Longitude"][1], point["Latitude"][1]))
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

    Node(
        Point(
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

    Edge(
        edge_dict["source"],
        edge_dict["target"],
        id,
        points
    )
end

function read_graph(file)
    nodes = []
    edges = []
    open(file) do f
        dict = parseGMLDict(f; log=true)["graph"]
        for node in dict["node"]
            push!(nodes, create_node(node))
        end
        for (i, edge) in enumerate(dict["edge"])
            push!(edges, create_edge(edge, i))
        end
    end

    nodes, edges
end


lat(node::Node) = node.point.lat
lon(node::Node) = node.point.lon

lat(point::Point) = point.lat
lon(point::Point) = point.lon


function vertex_str(n::Node; label=true, id=false)

    @assert !(label && id) "Label and id should not be used together"

    r_label = label ? n.label : (id ? n.id : "")
    r_label = replace(r_label, '_' => "\\_")
    "\\Vertex[y=$(n.point.lat), x=$(n.point.lon), label=$(r_label)]{$(n.id)}"
end

function node_str(p::Point, label; draw=false)
    "\\Vertex[y=$(p.lat), x=$(p.lon), shape = coordinate]{$(label)}"
    #"\\node[x=$(p.lat), y=$(p.lon), $(draw ? "" : "draw = none")]($(label)){};"
end

function get_point(p::Point)
    "{$(p.lon), $(p.lat)}"
end

function get_points(ps::Vector{Point})
    join(get_point.(ps), ',')
end

function edge_str(e::Edge; id=false, bend=0, style="solid", color="black")::String
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

nodes, edges = read_graph("graphs/italy6/italy6.gml")

lines = graph_to_tikz_net(nodes, edges)
for l in lines
    println(l)
end
