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


end