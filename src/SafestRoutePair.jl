module SafestRoutePair

include("graph_utils.jl")

using .GraphUtils
import .GraphUtils: weighted_graph_from_mat, calc_edges_from_nodes, create_indep_graph

include("readers/reader.jl")
include("algorithms/aco.jl")
include("algorithms/ga.jl")
include("algorithms/naive_complex.jl")
include("plotters/tikz_plotter.jl")


using Graphs
using Folds
using .ACO
import .ACO: ACOSettings, AcoRunSettings, ACO_preprocessing, calc_fitness

using .GA
import .GA: GeneticSettings, GaRunSettings, genetic, mutate_permute, mutate_random, npoint_crossover_naive, one_point_crossover_naive, crossover_roulette, calc_fitness_paths, calc_fitness_sets

using .NaiveComplex
import .NaiveComplex: naive_complex, NaiveGreedySettings

using .GraphReader
import .GraphReader: read_graph_and_failure, read_graph_with_positions, GeoPoint, NodeWithGeoPoint, EdgeWithGeoPoints, lat, lon

using .TIKZPlotter
import .TIKZPlotter: graph_to_tikz_net



export ACO_preprocessing, calc_fitness, ACOSettings, AcoRunSettings
export GeneticSettings, GaRunSettings, genetic, mutate_permute, mutate_random, npoint_crossover_naive, one_point_crossover_naive, crossover_roulette, calc_fitness_paths, calc_fitness_sets
export read_graph_and_failure, read_graph_with_positions
export graph_to_tikz_net, lat, lon



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


function find_id_of_single_edge(edge, fp_edges)
    for (i, edges) in enumerate(fp_edges)
        if length(edges) == 1 && (edge == edges[1] || (snd(edge), fst(edge)) == edges[1])
            return i
        end
    end

    return -1
end

function is_edge_in_list(list, edge)
    for e in list
        if (fst(e) == fst(edge) && snd(e) == snd(edge)) || (fst(e) == snd(edge) && snd(e) == fst(edge))
            return true
        end
    end
    false
end

function get_set_of_edge(p1_edges, p2_edges, edges)
    if length(edges) > 1
        return 0
    end
    any_in_p1 = any([is_edge_in_list(p1_edges, edge) for edge in edges])
    any_in_p2 = any([is_edge_in_list(p2_edges, edge) for edge in edges])

    if any_in_p1 && any_in_p2
        0
    elseif any_in_p1
        2
    elseif any_in_p2
        1
    else
        rand([1, 2])
    end
end

function result_to_ga_chromosome((p1, p2), fp_edges)
    p1_edges = collect(zip(p1, p1[2:end]))
    p2_edges = collect(zip(p2, p2[2:end]))

    [get_set_of_edge(p1_edges, p2_edges, edges) for edges in fp_edges]
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
    @time "ACO: $(from), $(to)" res1, _val1 = SafestRoutePair.ACO_preprocessing(acoS, runS1, gcfp.g, gcfp.cfps, gcfp.cfp_edges; logging_file=logging_file, use_folds=use_folds)

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
prob_limit ::Real - (Optional) FP values below this value are ignored in order to reduce chromosome size (0 means disabled)
edge_limit ::Int - (Optional) FP-s containing more edges than this value are ignored in order to reduce chromosome size (0 means disabled), if prob_limit is provided, this will be ignored
fp_depend ::Vector{Vector{Integer}} - (Optional) dependencies between FP-s, if an FP set contains two or more edges, we consider it to be dependent on all FP-s which only contain those edges. This only works with calc_fitness_sets. 
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
function safest_route_pair_ga(gcfp::GraphWithFPandCFP, from, to, gaS::GeneticSettings; logging_file="", use_folds=true, fp_depend=[])

    runS1 = []
    if gaS.prob_limit == 0 && gaS.edge_limit == 0
        runS1 = SafestRoutePair.GaRunSettings(gcfp.g, gcfp.fps, gcfp.fp_edges, fp_depend, gcfp.cfps, gcfp.cfp_edges, from, to)
    elseif gaS.prob_limit != 0
        # 0.0001
        indices = gcfp.fps .>= gaS.prob_limit
        new_fps = deepcopy(gcfp.fps[indices])
        new_fp_edges = deepcopy(gcfp.fp_edges[indices])
        runS1 = GaRunSettings(gcfp.g, new_fps, new_fp_edges, gcfp.cfps, gcfp.cfp_edges, from, to)
    else
        indices = length.(gcfp.fp_edges) .<= gaS.edge_limit
        new_fps = deepcopy(gcfp.fps[indices])
        new_fp_edges = deepcopy(gcfp.fp_edges[indices])
        runS1 = GaRunSettings(gcfp.g, new_fps, new_fp_edges, gcfp.cfps, gcfp.cfp_edges, from, to)
    end



    chromosomes = [[length(es) == 1 ? rand([0, 1, 2]) : 0 for es in runS1.fp_edges] for _ in 1:gaS.populationSize] #new_fp_edges

    if gaS.seed_naive
        naive_result, _ = safest_route_pair_naive(gcfp, from, to)
        naive_chromosome = result_to_ga_chromosome(naive_result, runS1.fp_edges)
        for i in 1:Int(floor(gaS.populationSize / 10))
            chromosomes[i] = deepcopy(naive_chromosome)
        end

    end

    res1 = SafestRoutePair.genetic(runS1, gaS, chromosomes; logging_file=logging_file, use_folds=use_folds)


    #@show res1
    val1 = calc_availability(res1, gcfp.fps, gcfp.fp_edges)
    #@show val1

    (res1, val1)
end

"""

Runs the naive algorithm to find the safest route pair between `from` and `to` in graph `gcfp`.

Parameters (most types are not specified, are only listed for reference)

```
gcfp  ::GraphWithFPandCFP - graph and probabilities
from  ::Int64 - the index of the starting node (julia indexes from 1)
to    ::Int64 - the index of the target node (julia indexes from 1)
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
function safest_route_pair_naive(gcfp::GraphWithFPandCFP, from, to)

    @time "Naive $(from), $(to)" res1 = SafestRoutePair.naive_complex(adjacency_matrix(gcfp.g), NaiveGreedySettings(from, to, gcfp.cfps, gcfp.cfp_edges))


    #@show res1
    val1 = calc_availability(res1, gcfp.fps, gcfp.fp_edges)
    #@show val1

    (res1, val1)
end



# Extracted in order to prevent unnecessary code duplication
function get_node_pairs(g, undirected::Bool, limit_pairs=0)
    if undirected
        if limit_pairs != 0
            return collect(Iterators.flatten([[(i, j) for j = (i+1):nv(g)] for i = 1:nv(g)]))[1:limit_pairs]
        end
        return collect(Iterators.flatten([[(i, j) for j = (i+1):nv(g)] for i = 1:nv(g)]))
    end

    if limit_pairs != 0
        return collect(filter(((x, y),) -> x != y, Iterators.flatten([[(i, j) for j = 1:nv(g)] for i = 1:nv(g)])))[1:limit_pairs]
    end

    collect(filter(((x, y),) -> x != y, Iterators.flatten([[(i, j) for j = 1:nv(g)] for i = 1:nv(g)])))
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
limit_pairs :: Integer - (Optional) how many pairs should be tested (0 means all pairs)
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
    acoS::ACOSettings=ACOSettings(10, 0.5, 100, 0.3, 0.1, 200, 10),
    logging_file="",
    use_folds=true,
    undirected=true,
    limit_pairs=0,
)

    node_pairs = get_node_pairs(gcfp.g, undirected, limit_pairs)


    run_algorithm(
        ((x, y),) -> ((x, y), safest_route_pair_aco(gcfp, x, y, acoS; logging_file=(logging_file != "" ? "$(logging_file)_$(x)_$(y).csv" : ""), use_folds=use_folds)),
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
limit_pairs :: Integer - (Optional) how many pairs should be tested (0 means all pairs)
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
    gaS::GeneticSettings=GeneticSettings(25, 0.1, 0.9, 0.5, npoint_crossover_naive, mutate_permute, 100, calc_fitness_sets, 0, 0, false),
    logging_file="",
    use_folds=false,
    undirected=true,
    limit_pairs=0,
)

    node_pairs = get_node_pairs(gcfp.g, undirected, limit_pairs)


    fp_depend = [[find_id_of_single_edge(edge, gcfp.fp_edges) for edge in edges] for edges in gcfp.fp_edges]

    run_algorithm(
        ((x, y),) -> ((x, y), safest_route_pair_ga(gcfp, x, y, gaS; logging_file=(logging_file != "" ? "$(logging_file)_$(x)_$(y).csv" : ""), use_folds=use_folds, fp_depend)),
        node_pairs,
        use_folds
    )
end

"""
Runs the genetic algorithm for all node pairs (by default undirected only, but configurable), where the nodes are distinct.

Parameters (most types are not specified, are only listed for reference)
```
gcfp  ::GraphWithFPandCFP - graph and probabilities
use_folds ::Bool - (Optional) should the algorithm use the Folds package to parallelize some parts (by default: true)
undirected :: Bool - (Optional) should the algorithm only look in one direction when searching a path between two nodes (by default: true)
limit_pairs :: Integer - (Optional) how many pairs should be tested (0 means all pairs)
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
function safest_route_pairs_all_naive(
    gcfp::GraphWithFPandCFP;
    use_folds=false,
    undirected=true,
    limit_pairs=0,
)

    node_pairs = get_node_pairs(gcfp.g, undirected, limit_pairs)


    run_algorithm(
        ((x, y),) -> ((x, y), safest_route_pair_naive(gcfp, x, y)),
        node_pairs,
        use_folds
    )
end

## gaS::GeneticSettings=GeneticSettings(50, 0.1, 0.9, 0.5, crossover_roulette, mutate, 200),

end # module AntSafestRoutePair
