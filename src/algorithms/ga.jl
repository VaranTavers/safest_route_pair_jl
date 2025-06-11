module GA

using Graphs
using SimpleWeightedGraphs
using Folds
using CSV
using DataFrames
using Base.Iterators
using ..GraphUtils

import .GraphUtils: weighted_graph_from_mat, calc_edges_from_nodes, create_indep_graph

#=
prob_limit ::Real - (Optional) FP values below this value are ignored in order to reduce chromosome size (0 means disabled)
edge_limit ::Int - (Optional) FP-s containing more edges than this value are ignored in order to reduce chromosome size (0 means disabled), if prob_limit is provided, this will be ignored
calculate_dependencies:: Bool - (optional) should it calculate dependencies between FP-s, if an FP set contains two or more edges, we consider it to be dependent on all FP-s which only contain those edges. This only works with calc_fitness_sets.
=#
struct GeneticSettings
    populationSize::Any
    mutationRate::Any
    crossoverRate::Any
    elitRate::Any
    crossoverAlg::Any
    mutationAlg::Any
    numberOfIterations::Any
    fitnessCalc::Any
    prob_limit::Any
    edge_limit::Any
    seed_naive::Bool
end

struct GaRunSettings
    g::Any
    fps::Vector{Real}
    fp_edges::Vector{Vector{Tuple{Integer,Integer}}}
    fp_depend::Vector{Vector{Integer}}
    cfps::Vector{Real}
    cfp_edges::Vector{Vector{Tuple{Integer,Integer}}}
    source
    target
end

fst((x, _)) = x
snd((_, y)) = y


sample(weights) = findfirst(cumsum(weights) .> rand())


function subgraphs_by_partition(g, partition, fp_edges)
    g1 = deepcopy(g)
    for (p, edge_list) in zip(partition, fp_edges)
        if p == 1
            for (src, dst) in edge_list
                rem_edge!(g1, src, dst)
            end
        end
    end

    g2 = deepcopy(g)
    for (p, edge_list) in zip(partition, fp_edges)
        if p == 2
            for (src, dst) in edge_list
                rem_edge!(g2, src, dst)
            end
        end
    end

    g1, g2
end

function vec_from_parents(end_p, parents, dist)
    if dist == Inf
        return []
    end
    res = zeros(Int64, length(parents))
    res[1] = end_p
    i = 1
    while (parents[res[i]] != 0)
        res[i+1] = parents[res[i]]
        i += 1
    end

    reverse(res[res.!=0])
end

function partition_to_paths(g_indep, partition, fp_edges, source, target)
    g1, g2 = subgraphs_by_partition(g_indep, partition, fp_edges)

    D1 = dijkstra_shortest_paths(g1, source)
    D2 = dijkstra_shortest_paths(g2, source)

    vec_from_parents(target, D1.parents, D1.dists[target]), vec_from_parents(target, D2.parents, D2.dists[target])
end

# Basic mutation that permutates a value in the 0->1->2->0 order on one of the chromosomes
function mutate_permute(partition)
    to_mutate = rand(1:length(partition))
    partition2 = deepcopy(partition)
    partition2[to_mutate] = (partition[to_mutate] + 1) % 3

    partition2
end

function mutate_random(partition)
    to_mutate = rand(1:length(partition))
    partition2 = deepcopy(partition)
    partition2[to_mutate] = rand(0:2)

    partition2
end

function mutate_random_multiple(partition, gene_mut_prob=0.01)
    [rand() < gene_mut_prob ? rand(0:2) : x for x in partition]
end

# Basic crossover that chooses each cromosome randomly from one of the parents
function npoint_crossover_naive(partition1, partition2)
    [rand((1, 2)) == 1 ? partition1[i] : partition2[i] for i in eachindex(partition1)]
end

function one_point_crossover_naive(partition1, partition2)
    point = rand(1:length(partition1))

    res = deepcopy(partition1[1:point])

    append!(res, partition2[point+1:end])

    res
end

function crossover_roulette(_g, chromosomes, fitness, crossoverAlg)
    rouletteWheel = ones(length(chromosomes))

    if maximum(minimum.(fitness)) > 0
        rouletteWheel = minimum.(fitness) ./ sum(minimum.(fitness))
    end
    crossoverAlg(
        chromosomes[sample(rouletteWheel)],
        chromosomes[sample(rouletteWheel)]
    )
end

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


function calc_fitness_paths(solution, g_indep, runS::GaRunSettings) # _complex

    paths = partition_to_paths(g_indep, solution, runS.fp_edges, runS.source, runS.target)


    if fst(paths) == [] && snd(paths) == []
        return -20
    elseif fst(paths) == []
        return -10
    elseif snd(paths) == []
        return -10
    end

    calc_availability(paths, runS.fps, runS.fp_edges)
end




# Alternative fitness function: sum of FP-s from both sets should be maximal
function calc_fitness_sets(solution, g_indep, runS::GaRunSettings)  #_simple

    paths = partition_to_paths(g_indep, solution, runS.fp_edges, runS.source, runS.target)


    if fst(paths) == [] && snd(paths) == []
        return -20
    elseif fst(paths) == []
        return -10
    elseif snd(paths) == []
        return -10
    end

    solution2 = deepcopy(solution)

    if length(runS.fp_depend) === length(runS.fp_edges)
        edge_comps = runS.fp_depend

        for (i, (sol, comp)) in enumerate(zip(solution, edge_comps))
            l = length(comp)
            if l > 1 && sol == 0
                vals = [solution[j] for j in comp]
                if all(vals .== 1)
                    solution2[i] = 1
                elseif all(vals .== 2)
                    solution2[i] = 2
                end
            end
        end
    end

    sum(runS.fps[solution2.!=0])
end



function genetic(runS::GaRunSettings, gaS::GeneticSettings, chromosomes; logging_file="", use_folds=true, logging_depth="iteration_best")

    # Initializing values and functions for later use
    g_indep = create_indep_graph(runS.g, runS.cfps, runS.cfp_edges)

    calcFitness(x) = gaS.fitnessCalc(x, g_indep, runS)
    runMutation(x) = rand() < gaS.mutationRate ? gaS.mutationAlg(x) : x
    chromosomes::Vector{Vector{Int64}} = deepcopy(chromosomes)
    n_c = length(chromosomes)

    # Initializing global maximum as one of the given chromosome
    maxVal = calcFitness(chromosomes[1])
    maxVec = copy(chromosomes[1])

    fitness = []
    if use_folds
        fitness = Folds.map(calcFitness, chromosomes)
    else
        fitness = collect(map(calcFitness, chromosomes))
    end

    #@show runS
    for i = 1:gaS.numberOfIterations
        #@show maxVal, maxVec
        #@show chromosomes
        if maxVal == Inf
            @show "GA: Early exit, due to unbreakable route"
            return partition_to_paths(g_indep, maxVec, runS.fp_edges, runS.source, runS.target)
        end
        if maximum(fitness) == Inf
            @show "GA: Early exit, due to unbreakable route"
            return partition_to_paths(g_indep, chromosomes[argmax(fitness)], runS.fp_edges, runS.source, runS.target)
        end
        # Creating p_c% new individuals with the crossover
        # operator, choosing parents based on fitness.
        newChromosomes = [
            crossover_roulette(runS.g, chromosomes, fitness, gaS.crossoverAlg) for
            _ = 1:Int(ceil(n_c * gaS.crossoverRate))
        ]

        newFitness = [0.0 for _ = 1:length(newChromosomes)]

        # Add them to the chromosome pool
        append!(chromosomes, newChromosomes)
        append!(fitness, newFitness)


        # Mutating individuals
        chromosomes = collect(map(x -> runMutation(x), chromosomes))

        # Recalculating fitness for new individuals
        if use_folds
            fitness = Folds.map(calcFitness, chromosomes)
        else
            fitness = collect(map(calcFitness, chromosomes))
        end
        # Sorting fitness scores
        fitnessSorted = sortperm(fitness, rev=true, by=minimum)

        fitnessMaxVal = deepcopy(fitness[fitnessSorted[1]])
        fitnessMaxVec = deepcopy(chromosomes[fitnessSorted[1]])

        # Choosing the elit
        elitNumber = Int(ceil(gaS.populationSize * gaS.elitRate))
        elitChromosomes = deepcopy(chromosomes[fitnessSorted[1:elitNumber]])
        elitFitness = copy(fitness[fitnessSorted[1:elitNumber]])


        # Choosing the rest randomly from the others
        restNumber = gaS.populationSize - elitNumber
        restIds = [rand(fitnessSorted[elitNumber+1:end]) for _ = 1:restNumber]
        restChromosomes = map(x -> copy(chromosomes[x]), restIds)
        restFitness = map(x -> fitness[x], restIds)

        chromosomes = vcat(elitChromosomes, restChromosomes)
        fitness = vcat(elitFitness, restFitness)

        if minimum(fitnessMaxVal) > minimum(maxVal)
            maxVec = deepcopy(fitnessMaxVec)
            maxVal = deepcopy(fitnessMaxVal)
        end

        if logging_file != ""
            logdf = []
            if logging_depth == "iteration_best"
                logdf = DataFrames.DataFrame(i=i, maxVec=[maxVec === nothing ? [] : maxVec], maxVal=[maxVal === nothing ? -1 : maxVal])
            elseif logging_depth == "all"
                logdf = DataFrames.DataFrame(i=i, vecs=["$(chromosomes)"], vals=["$(fitness)"])
            end
            CSV.write(logging_file, logdf; append=true)
        end
    end
    #@show maxVal
    #@show maxVec


    p1, p2 = partition_to_paths(g_indep, maxVec, runS.fp_edges, runS.source, runS.target)

    if p1 == []
        p1 = p2
    end
    if p2 == []
        p2 = p1
    end
    if p1 == [] && p2 == []
        @show "GA: No route found?"
    end
    (p1, p2)
end

end
