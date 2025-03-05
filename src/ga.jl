module GA

using Graphs
using SimpleWeightedGraphs
using Folds
using Base.Iterators

struct GeneticSettings
    populationSize::Any
    mutationRate::Any
    crossoverRate::Any
    elitRate::Any
    crossoverAlg::Any
    mutationAlg::Any
    numberOfIterations::Any
end

struct GaRunSettings
    g::Any
    fps::Vector{Real}
    fp_edges::Vector{Vector{Tuple{Integer,Integer}}}
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
function mutate(partition)
    to_mutate = rand(1:length(partition))
    partition[to_mutate] = (partition[to_mutate] + 1) % 3
end

# Basic crossover that chooses each cromosome randomly from one of the parents
function single_crossover_naive(partition1, partition2)
    [rand((1, 2)) == 1 ? partition1[i] : partition2[i] for i in eachindex(partition1)]
end

function crossover_roulette(_g, chromosomes, fitness)
    rouletteWheel = ones(length(chromosomes))

    if maximum(minimum.(fitness)) > 0
        rouletteWheel = minimum.(fitness) ./ sum(minimum.(fitness))
    end
    single_crossover_naive(
        chromosomes[sample(rouletteWheel)],
        chromosomes[sample(rouletteWheel)]
    )
end

function calc_edges_from_nodes(x)
    collect(zip(x, x[2:end]))
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

    @show path_a
    @show path_b

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

function calc_fitness(solution, g_indep, runS::GaRunSettings)

    paths = partition_to_paths(g_indep, solution, runS.fp_edges, runS.source, runS.target)

    @show paths

    if fst(paths) == [] && snd(paths) == []
        return -2000
    elseif fst(paths) == []
        return -1000
    elseif snd(paths) == []
        return -1000
    end

    calc_availability(paths, runS.fps, runS.fp_edges)
end

function weighted_graph_from_mat(mat)
    n, _ = size(mat)
    g = SimpleWeightedGraph(n)

    for i = 1:n
        for j = 1:n
            add_edge!(g, i, j, mat[i, j])
        end
    end

    g
end


# Something wrong here
function create_indep_graph(g, cfps::Vector{Real}, cfp_edges::Vector{Vector{Tuple{Integer,Integer}}})
    edges_mat = adjacency_matrix(g)
    hits = ones(size(edges_mat)) .* 0.0000001

    # Adds single edge CFP-s as weights to the graph as weights
    for (break_val, edges) in zip(cfps, cfp_edges)
        if length(edges) == 1
            s, e = edges[1]
            hits[s, e] = break_val
            hits[e, s] = break_val
        end
    end

    weighted_graph_from_mat(hits .* edges_mat)
end

function genetic(runS::GaRunSettings, gaS::GeneticSettings, chromosomes; logging_file="", use_folds=true)
    # Initializing values and functions for later use
    g_indep = create_indep_graph(runS.g, runS.cfps, runS.cfp_edges)

    calcFitness(x) = calc_fitness(x, g_indep, runS)
    runMutation(x) = rand() < gaS.mutationRate ? gaS.mutationAlg(x) : x
    chromosomes = deepcopy(chromosomes)
    n_c = length(chromosomes)

    # Initializing global maximum as one of the given chromosome
    maxVal = calcFitness(chromosomes[1])
    maxVec = copy(chromosomes[1])

    # Initializing logging
    logs = []

    fitness = []
    if use_folds
        fitness = Folds.map(calcFitness, chromosomes)
    else
        fitness = collect(map(calcFitness, chromosomes))
    end
    @show fitness

    for i = 1:gaS.numberOfIterations

        if maxVal == Inf
            @show "GA: Early exit, due to unbreakable route"
            return partition_to_paths(g_indep, maxVec, runS.fp_edges, runS.source, runS.target), logs
        end
        if maximum(fitness) == Inf
            @show "GA: Early exit, due to unbreakable route"
            return partition_to_paths(g_indep, chromosomes[argmax(fitness)], runS.fp_edges, runS.source, runS.target), logs
        end
        # Creating p_c% new individuals with the crossover
        # operator, choosing parents based on fitness.
        newChromosomes = [
            gaS.crossoverAlg(runS.g, chromosomes, fitness) for
            _ = 1:Int(ceil(n_c * gaS.crossoverRate))
        ]

        newFitness = [0.0 for _ = 1:length(newChromosomes)]

        @show fitness
        @show newFitness
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
            logRow = (i, maxVal, maxVec)
            push!(logs, logRow)
        end
    end

    partition_to_paths(g_indep, maxVec, runS.fp_edges, runS.source, runS.target), logs
end

end