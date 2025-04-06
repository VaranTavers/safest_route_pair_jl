module GraphReader
using XML
using Graphs
using SimpleWeightedGraphs
using GraphIO

include("gml_reader.jl")
using .GMLReader

import .GMLReader: GeoPoint, NodeWithGeoPoint, EdgeWithGeoPoints, read_graph_with_positions, lat, lon
export GeoPoint, NodeWithGeoPoint, EdgeWithGeoPoints, read_graph_with_positions, lat, lon

function getTagFromChildren(node, tagName)
    for x in children(node)
        if tag(x) == tagName
            return x
        end
    end

    nothing
end

function parseEdgeRow(string::AbstractString)
    vals = split(string, '(')
    src_dst = split(vals[2], ',')
    if contains(src_dst[1], ':')
        src = split(src_dst[1], ':')
        dst = split(src_dst[2], ':')
    else
        src = [src_dst[1]]
        dst = split(src_dst[2], ')')
    end

    (parse(Int64, src[1]) + 1, parse(Int64, dst[1]) + 1)
end

function parseEdgeNumbers(string::AbstractString)
    delim = "\n"
    if contains(string, "\r")
        delim = "\r\n"
    end
    rows = split(string, delim)

    collect(map(parseEdgeRow, rows))
end

function parseFailureState(node::Node)
    prob = parse(Float64, simplevalue(getTagFromChildren(node, "Probability")))
    if !is_simple(getTagFromChildren(node, "Edges"))
        return (0, [])
    end
    edges_string = simplevalue(getTagFromChildren(node, "Edges"))
    #@show edges_string

    edge_values = parseEdgeNumbers(edges_string)
    prob, edge_values
end

function parseFailureFile(failure_file)
    doc = read(failure_file, Node)

    tree = doc[end]

    if tree.tag == "Failure_State_Distribution"
        return [parseFailureState(x) for x in children(tree)[2:end]]
    elseif tree.tag == "simulation"
        return [parseFailureState(x) for x in children(children(tree)[6])]
    end

end

#=
    probs = []

    for x in children(tree)
        if x.tag == "Failure_State"

        else if x.tag == "PSRLGList"
            return 
        end
    end
=#



function get_first_graph_name(graph_file)
    ret = "unknown"
    open(graph_file) do f
        while !eof(f)
            s = readline(f)
            if contains(s, "label")
                parts = collect(split(s, "\""))
                ret = parts[2]
                break
            end
        end

    end

    "$(ret)"
end


function read_graph_and_failure(graph_file, cfp_file, fp_file)

    cfp_s = parseFailureFile(cfp_file)
    fp_s = parseFailureFile(fp_file)



    graph_name = get_first_graph_name(graph_file)
    # "../interroute_j.gml", "interroute_v2"
    g = loadgraph(graph_file, graph_name, CustomGMLFormat())


    #@show collect(src.(edges(g)))
    #@show collect(dst.(edges(g)))

    sample(weights) = findfirst(cumsum(weights) .> rand())

    fst((x, _)) = x
    snd((_, y)) = y

    #@show sum(fst.(probs))
    #@show sample(fst.(probs))

    #@show fst.(probs), snd.(probs)

    g, fst.(cfp_s), snd.(cfp_s), fst.(fp_s), snd.(fp_s)

end
end