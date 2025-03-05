module SafestRoutePair

include("reader.jl")
include("aco.jl")
include("ga.jl")

using Graphs
using Folds
using .ACO
import .ACO: ACOSettings, AcoRunSettings, ACO_preprocessing, calc_fitness

using .GA
import .GA: GeneticSettings, GaRunSettings, genetic, mutate, crossover_roulette

using .GraphReader
import .GraphReader: read_graph_and_failure


export read_graph_and_failure
export ACO_preprocessing, calc_fitness
export ACOSettings, AcoRunSettings

export GeneticSettings, GaRunSettings, genetic, mutate, crossover_roulette


"""
Type that packs the graph together with the breakage probabilities of certain edge sets.

Properties:
```
    g ::AbstractGraph - the graph on which we run the code
    fps ::Vector{Real} - the probabilities of breakage (FP) of the associated edge group
    fp_edges ::Vector{Vector{Tuple{Integer,Integer}}} - the edge groups that have an associated FP
    cfps ::Vector{Real} - the probabilities of breakage (CFP) of the associated edge group
    cfp_edges ::Vector{Vector{Tuple{Integer,Integer}}} - the edge groups that have an associated CFP
```
"""

struct GraphWithFPandCFP
    g
    fps::Vector{Real}
    fp_edges::Vector{Vector{Tuple{Integer,Integer}}}
    cfps::Vector{Real}
    cfp_edges::Vector{Vector{Tuple{Integer,Integer}}}
end


function calc_edges_from_nodes(x)
    collect(zip(x, x[2:end]))
end

fst((x, _)) = x
snd((_, y)) = y

# Checks if any y-s are present in x-s
function path_intersects(xs, ys)
    for y in ys

        if findfirst(fst.(xs) .== fst(y) .&& snd.(xs) .== snd(y)) !== nothing ||
           findfirst(fst.(xs) .== snd(y) .&& snd.(xs) .== fst(y)) !== nothing
            return true
        end
    end

    false
end

"""
Calculates the availability of the path pair given in `path1` and `path2` using the FP values in `fps` and the edge sets to which they are associated.
"""
function calc_availability((path1, path2), fps::Vector{Real}, fp_edges::Vector{Vector{Tuple{Integer,Integer}}})
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


"""

Runs the ACO algorithm to find the safest route pair between `from` and `to` in graph `gcfp`.

Parameters (most types are not specified, are only listed for reference)

```
gcfp  ::GraphWithFPandCFP - graph and probabilities
from  ::Int64 - the index of the starting node (julia indexes from 1)
to    ::Int64 - the index of the target node (julia indexes from 1)
acoS  ::ACOSetting - the parameters for the ACO algorithm
logging_file ::String - (Optional) name of the file where we save the logs, leave empty for no logging
use_folds ::Bool - (Optional) should the algorithm use the Folds package to parallelize some parts
```

Return value:
```
Tuple{Vector{Int64}, Float64}
(
  res ::Tuple{Vector{Int64}, Vector{Int64}} - two list of nodes describing the two paths
  val ::Float64 - the negative logarithm of the breakage probability (to convert it back do: 1 - exp(-val))
)
```
"""
function safest_route_pair_aco(gcfp::GraphWithFPandCFP, from, to, acoS; logging_file="", use_folds=true)
    #@show (from, to)
    #@show Matrix(adjacency_matrix(g))



    runS1 = SafestRoutePair.AcoRunSettings(from, to, from + nv(gcfp.g), gcfp.fps, gcfp.fp_edges)
    res1, _val1 = SafestRoutePair.ACO_preprocessing(acoS, runS1, gcfp.g, gcfp.cfps, gcfp.cfp_edges; logging_file=logging_file, use_folds=use_folds)

    val1 = calc_availability(res1, gcfp.fps, gcfp.fp_edges)

    (res1, val1)
end

"""

Runs the genetic algorithm to find the safest route pair between `from` and `to` in graph `gcfp`.

Parameters (most types are not specified, are only listed for reference)

```
gcfp  ::GraphWithFPandCFP - graph and probabilities
from  ::Int64 - the index of the starting node (julia indexes from 1)
to    ::Int64 - the index of the target node (julia indexes from 1)
gaS   ::GeneticSettings - the parameters for the genetic algorithm
logging_file ::String - (Optional) name of the file where we save the logs, leave empty for no logging
use_folds ::Bool - (Optional) should the algorithm use the Folds package to parallelize some parts
```

Return value:
```
Tuple{Vector{Int64}, Float64}
(
  res ::Tuple{Vector{Int64}, Vector{Int64}} - two list of nodes describing the two paths
  val ::Float64 - the negative logarithm of the breakage probability (to convert it back do: 1 - exp(-val))
)
```
"""
function safest_route_pair_ga(gcfp::GraphWithFPandCFP, from, to, gaS::GeneticSettings; logging_file="", use_folds=true)

    runS1 = SafestRoutePair.GaRunSettings(gcfp.g, gcfp.fps, gcfp.fp_edges, gcfp.cfps, gcfp.cfp_edges, from, to)
    chromosomes = [[rand([0, 1, 2]) for _ in 1:length(gcfp.fps)] for _ in 1:gaS.populationSize]
    @show from, to
    res1, _logs = SafestRoutePair.genetic(runS1, gaS, chromosomes; logging_file=logging_file, use_folds=use_folds)

    val1 = calc_availability(res1, gcfp.fps, gcfp.fp_edges)

    (res1, val1)
end



# Extracted in order to prevent unnecessary code duplication
function get_node_pairs(g, undirected::Bool)
    if undirected
        return collect(Iterators.flatten([[(i, j) for j = (i+1):nv(g)] for i = 1:nv(g)]))
    end

    collect(filter(((x, y),) -> x != y, Iterators.flatten([[(i, j) for j = 1:nv(g)] for i = 1:nv(g)])))[1:1]
end

# Extracted in order to prevent unnecessary code duplication
function run_algorithm(f, node_pairs, use_folds::Bool)
    if use_folds
        return Folds.map(f, node_pairs)
    end

    collect(map(f, node_pairs))

end

"""
Runs the ACO algorithm for all node pairs (by default undirected only, but configurable), where the nodes are distinct.

Parameters (most types are not specified, are only listed for reference)
```
gcfp  ::GraphWithFPandCFP - graph and probabilities
acoS  ::ACOSetting (Optional) - the parameters for the ACO algorithm
logging_file ::String - (Optional) name of the file where we save the logs, leave empty for no logging (by default: "")
use_folds ::Bool - (Optional) should the algorithm use the Folds package to parallelize some parts (by default: true)
undirected :: Bool - (Optional) should the algorithm only look in one direction when searching a path between two nodes (by default: true)
```

Return value:

Pairs of routes and the negative logarithms of breakage probabilities.
```
Vector{Tuple{(Vector{Int64}, Float64})}}
[(
  res ::Tuple{Vector{Int64}, Vector{Int64}} - two list of nodes describing the two paths
  val ::Float64 - the negative logarithm of the breakage probability (to convert it back do: 1 - exp(-val))
)]
```
"""
function safest_route_pairs_all_aco(
    gcfp::GraphWithFPandCFP;
    acoS::ACOSettings=ACOSettings(1, 0.5, 100, 0.3, 0.1, 200, 1),
    logging_file="",
    use_folds=true,
    undirected=true,
)

    node_pairs = get_node_pairs(gcfp.g, undirected)


    run_algorithm(
        ((x, y),) -> ((x, y), safest_route_pair_aco(gcfp, x, y, acoS; logging_file=(logging_file != "" ? "$(logging_file)_$(x)_$(y).csv" : ""))),
        node_pairs,
        use_folds
    )
end

"""
Runs the genetic algorithm for all node pairs (by default undirected only, but configurable), where the nodes are distinct.

Parameters (most types are not specified, are only listed for reference)
```
gcfp  ::GraphWithFPandCFP - graph and probabilities
gaS  ::GeneticSettings (Optional) - the parameters for the ACO algorithm
logging_file ::String - (Optional) name of the file where we save the logs, leave empty for no logging (by default: "")
use_folds ::Bool - (Optional) should the algorithm use the Folds package to parallelize some parts (by default: true)
undirected :: Bool - (Optional) should the algorithm only look in one direction when searching a path between two nodes (by default: true)
```

Return value:

Pairs of routes and the negative logarithms of breakage probabilities.
```
Vector{Tuple{(Vector{Int64}, Float64})}}
[(
  res ::Tuple{Vector{Int64}, Vector{Int64}} - two list of nodes describing the two paths
  val ::Float64 - the negative logarithm of the breakage probability (to convert it back do: 1 - exp(-val))
)]
```
"""
function safest_route_pairs_all_ga(
    gcfp::GraphWithFPandCFP;
    gaS::GeneticSettings=GeneticSettings(10, 0.1, 0.9, 0.5, crossover_roulette, mutate, 1),
    logging_file="",
    use_folds=false,
    undirected=true,
)

    node_pairs = get_node_pairs(gcfp.g, undirected)


    run_algorithm(
        ((x, y),) -> ((x, y), safest_route_pair_ga(gcfp, x, y, gaS; logging_file=(logging_file != "" ? "$(logging_file)_$(x)_$(y).csv" : ""))),
        node_pairs,
        use_folds
    )
end

## gaS::GeneticSettings=GeneticSettings(50, 0.1, 0.9, 0.5, crossover_roulette, mutate, 200),

end # module AntSafestRoutePair
