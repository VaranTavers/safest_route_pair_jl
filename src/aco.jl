module ACO
begin
    using Graphs
    using DataFrames
    using CSV
    using SimpleWeightedGraphs
    using Folds
    using ConcurrentCollections
end


struct ACOSettings
    α::Real
    β::Real
    number_of_ants::Integer
    ρ::Real
    ϵ::Real
    max_number_of_iterations::Integer
    starting_pheromone_ammount::Real
    ACOSettings(α, β, n_a, ρ, ϵ, max_i, start_ph) = new(α, β, n_a, ρ, ϵ, max_i, start_ph)
end

struct AcoRunSettings
    start_p::Integer
    mid_p::Integer
    end_p::Integer
    fps::Vector{Real}
    fp_edges::Vector{Vector{Tuple{Integer,Integer}}}
end


mutable struct ACOInner
    n::Any
    η::Any
    τ::Any
end

sample(weights) = findfirst(cumsum(weights) .> rand())



# Get chosen point
function get_chosen_point(pM, i, r)
    if maximum(pM[i, :]) == 0
        return i
    end

    findfirst(pM[i, :] .> r)
end



# Constructs a new solution in the form of a list of nodes visited: [1, 2, 3] 
function generate_s(n, pM, start_p, end_p)
    s = zeros(Int32, n)
    visited = [false for i ∈ 1:n]
    s[1] = start_p
    visited[start_p] = true
    i = 1
    while s[i] != end_p && i < n
        j = 0
        res = start_p
        while visited[res] && j < 100
            res = get_chosen_point(pM, s[i], rand())
            j += 1
        end
        # Stuck between a rock and a hard place
        if j == 100
            # @show "stuck"
            # @show visited
            # @show s
            return nothing
        end

        # Append to solution
        s[i+1] = res
        visited[res] = true
        i += 1
    end

    # Did not reach the goal
    if i == n && s[i] != end_p
        return nothing
    end

    s
end

# Transform solution from list of nodes into list of edges (node pairs)
function calc_edges_from_nodes(x)
    collect(zip(x, x[2:end]))
end

fst((x, _)) = x
snd((_, y)) = y

# Checks if all y-s are present in x-s
function is_present(xs, ys)
    for y in ys
        if findfirst(fst.(xs) .== fst(y) .&& snd.(xs) .== snd(y)) === nothing &&
           findfirst(fst.(xs) .== snd(y) .&& snd.(xs) .== fst(y)) === nothing
            return false
        end
    end

    true
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

# Splits single route of mirrored graph into two routes on the original graph

function split_edges(runS::AcoRunSettings, s_edges)
    n = runS.end_p - runS.start_p
    first_path = true
    path_a = []
    path_b = []

    for e in s_edges
        if first_path
            push!(path_a, e)
        elseif snd(e) != 0 && fst(e) != runS.mid_p
            push!(path_b, e .- n)
        end
        if snd(e) == runS.mid_p
            first_path = false
        end
    end

    (path_a, path_b)
end

# Splits lists of nodes of mirrored graph into two lists on the original graph

function split_final(runS::AcoRunSettings, s_edges)
    n = runS.end_p - runS.start_p
    first_path = true
    path_a = []
    path_b = []


    for e in s_edges
        if first_path
            push!(path_a, e)
        elseif e != 0
            push!(path_b, e - n)
        end
        if e == runS.mid_p
            first_path = false
            #push!(path_b, e)
        end
    end

    (path_a, reverse(path_b))
end


function calc_fitness(runS::AcoRunSettings, solution; fitnessCache=Dict(), logging=false)
    if solution === nothing
        return 0
    end
    if solution in keys(fitnessCache)
        return fitnessCache[solution]
    end
    # Transform solution from list of nodes into list of edges (node pairs)
    s_edges = calc_edges_from_nodes(solution)

    path_a, path_b = split_edges(runS, s_edges)

    only_path_a = 0
    only_path_b = 0
    both_paths = 0

    for (prob, vec) in zip(runS.fps, runS.fp_edges)
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


    fitnessCache[solution] = -log(res)


    fitnessCache[solution]
end


# η should be -log(edge failure), scoring should be sieve tech
function ACO_test(vars::ACOSettings, runS::AcoRunSettings, η, τ; logging_file="", use_folds=true)
    #Set parameters and initialize pheromone traits.
    n, _ = size(η)

    sgb = [i for i ∈ 1:n]
    sgb_val = -1000
    τ_max = vars.starting_pheromone_ammount
    τ_min = 0
    fitnessCache = ConcurrentDict{Vector{Integer},Float64}()


    #@show runS.start_p, runS.end_p
    # While termination condition not met
    for i ∈ 1:vars.max_number_of_iterations
        #@show i
        #@time 
        begin
            # Precomputing this
            η_d_sq = η .^ vars.β

            # Construct new solution S
            # Precomputing the probabilities results in a 2s time improvement.
            probM = τ .^ vars.α .* η_d_sq

            #probM = [τ[i, j] ^ vars.α * η[i, j] ^ vars.β for i in 1:n, j in 1:n]
            probM ./= sum(probM, dims=2)

            probM = cumsum(probM, dims=2)

            S = []
            if use_folds
                S = Folds.map(
                    x -> generate_s(n, probM, runS.start_p, runS.end_p),
                    1:vars.number_of_ants,
                )
            else
                S = collect(map(
                    x -> generate_s(n, probM, runS.start_p, runS.end_p),
                    1:vars.number_of_ants,
                ))
            end



            fitness = collect(map(x -> calc_fitness(runS, x; fitnessCache, logging=false), S))


            # Update iteration best
            maxIndex = argmax(fitness)
            (sib, sib_val) = (deepcopy(S[maxIndex]), fitness[maxIndex])
            if sib_val > sgb_val && sib !== nothing
                sgb_val = copy(sib_val)
                sgb = copy(sib)

                # Compute pheromone trail limits
                τ_max = sgb_val / (1 - vars.ρ)
                τ_min = vars.ϵ * τ_max
                if sgb_val == Inf
                    @show "ACO: Early exit, due to unbreakable route"
                    return split_final(runS, sgb), sgb_val
                end

            elseif sib === nothing
                @show "Not a single ant has reached to the target", sib_val, sgb_val
            end
            if !isempty(logging_file)
                #@show logging_file
                logdf = DataFrames.DataFrame(sib=[sib === nothing ? [] : sib], sib_val=[sib_val === nothing ? -1 : sib_val], sgb=[sgb === nothing ? [] : sgb], sgb_val=[sgb_val === nothing ? -1 : sgb_val])
                CSV.write(logging_file, logdf; append=true)
            end
            # Update pheromone trails
            # TODO: test with matrix sum
            τ .*= (1 - vars.ρ)
            if sib !== nothing
                for (s, e) in zip(sib, filter(x -> x > 0, sib[2:end]))
                    τ[s, e] += sib_val
                end
            end
            τ = min.(τ, τ_max)
            τ = max.(τ, τ_min)
        end
    end

    if sgb_val == -1000
        return ([], []), -1000
    end
    split_final(runS, sgb), sgb_val
end



function mirror_graph(g, probs, edge_groups, target)
    g2 = deepcopy(g)

    n = nv(g2)

    add_vertices!(g2, n)
    for e in edges(g)
        add_edge!(g2, e.src + n, e.dst + n)
    end

    add_edge!(g2, target, target + n)

    edges2 = deepcopy(edge_groups)
    append!(edges2, map(x -> map(y -> y .+ n, x), edge_groups))

    probs2 = deepcopy(probs)
    append!(probs2, probs)



    g2, probs2, edges2

end


function ACO_preprocessing(vars::ACOSettings, runS::AcoRunSettings, g, cfps, cfp_edges; logging_file="", use_folds=true)

    g2, cfps2, cfp_edges2 = mirror_graph(g, cfps, cfp_edges, runS.mid_p)
    inc_mat = adjacency_matrix(g2)

    neglogprobs = -1 .* log.(cfps2)
    # newRunS = RunSettings(runS.start_p, runS.end_p, neglogprobs, runS.fp_edges)
    τ = ones(size(inc_mat)) * vars.starting_pheromone_ammount
    η = deepcopy(inc_mat) .* -log(0.00000001)

    for (i, edges) in enumerate(cfp_edges2)
        if length(edges) == 1
            (node1, node2) = edges[1]
            η[node1, node2] = neglogprobs[i]
            η[node2, node1] = neglogprobs[i]
        end
    end


    ACO_test(vars, runS, η, τ; logging_file=logging_file, use_folds=use_folds)
end

end
