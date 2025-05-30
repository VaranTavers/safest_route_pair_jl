using Pkg

Pkg.activate(".")

using SafestRoutePair
using Graphs
using DataFrames
using CSV
using Statistics
using Dates
using Folds

fst((x, _)) = x
snd((_, y)) = y

function toString(x)
    "$(x)"
end

function create_dfs_from_results(results)
    local_route_df = DataFrame()
    local_result_df = DataFrame()

    local_route_df[!, :nodes] = fst.(results[1])
    local_result_df[!, :nodes] = fst.(results[1])

    for (i, result) in enumerate(results)
        local_route_df[!, :x] = toString.(fst.(snd.(result)))
        rename!(local_route_df, :x => "run$(i)")
        local_result_df[!, :x] = snd.(snd.(result))
        rename!(local_result_df, :x => "run$(i)")
    end

    tmpMat = Matrix(local_result_df[:, 2:(number_of_runs+1)])
    means = mean(tmpMat, dims=2)
    stds = std(tmpMat, dims=2)
    mins = minimum(tmpMat, dims=2)
    maxs = maximum(tmpMat, dims=2)


    local_result_df[!, :mean] = means[:]
    local_result_df[!, :std] = stds[:]
    local_result_df[!, :min] = mins[:]
    local_result_df[!, :max] = maxs[:]

    local_route_df, local_result_df
end

include("./test_config.jl")



configurations_ACO = [
    # conf_name,                        n_p, mut, cro, elit, crossoverAlg, mutationAlg,       meme,  log,  iter 
    # Baselines
    ("ACO_DEFAULT", ACOSettings(α, β, nr_ants, ρ, ϵ, nr_gen, s_pheromone)) for
    α in aco_param_tuning_α, β in aco_param_tuning_β, nr_ants in aco_param_tuning_nr_ants,
    ρ in aco_param_tuning_ρ, ϵ in aco_param_tuning_ϵ, nr_gen in aco_param_tuning_nr_gen,
    s_pheromone in aco_param_tuning_starting_pheromone
]

configurations_GA = [
    ("GA_SIMPLER_GENES_Paths2", GeneticSettings(popSize, mutRate, crossRate, elitRate, crossFunc, mutFunc, nrIter, calcFit, probLimit, edgeLimit, seedNaive)) for
    popSize in ga_param_tuning_n_p, mutRate in ga_param_tuning_mut_rate, crossRate in ga_param_tuning_cro_rate,
    elitRate in ga_param_tuning_elit, nrIter in ga_param_tuning_nr_gen, crossFunc in ga_param_tuning_crossover,
    mutFunc in ga_param_tuning_mutation, calcFit in ga_param_tuning_fitness, probLimit in ga_param_prob_limit, edgeLimit in ga_param_edge_limit, seedNaive in ga_param_seed_naive
]


println("Alg_tester_cli started on $(Dates.now())")

if !isdir("./graphs")
    println("The graphs directory is missing.")
    exit(1)
end

folders = readdir("./graphs")

# Creating the necessary folders

if !isdir("./logs")
    println("Creating ./logs directory")
    mkdir("./logs")
end
if !isdir("./results")
    println("Creating ./results directory")
    mkdir("./results")
end

date_of_start = Dates.today()

version_on_date = 1

while isdir("results/$(date_of_start)_$(version_on_date)")
    println("Folder results/$(date_of_start)_$(version_on_date) already exists, moving on to next number")
    global version_on_date += 1
end

res_folder_name = "$(date_of_start)_$(version_on_date)"
if !isdir("logs/$(res_folder_name)")
    mkdir("logs/$(res_folder_name)")
end
if !isdir("results/$(res_folder_name)")
    mkdir("results/$(res_folder_name)")
end

# Reading graphs

gcfps::Vector{Tuple{String,SafestRoutePair.GraphWithFPandCFP}} = []

for (i, folder) in enumerate(folders)
    files = readdir("graphs/$(folder)")

    graph_files = collect(filter(x -> x[(end-3):end] == ".gml", files))
    if length(graph_files) < 1
        println("No GML file was found in $(folder). Skipping")
        continue
    end
    if length(graph_files) > 1
        println("Multiple GML files were found in $(folder). The program will skip all but $(graph_files[1])")
    end

    graph_file = graph_files[1]
    cfp_file = "cfp.xml"
    fp_file = "fp.xml"

    if !isfile("graphs/$(folder)/$(cfp_file)")
        println("No CFP file found in $(folder). Skipping")
        continue
    end
    if !isfile("graphs/$(folder)/$(fp_file)")
        println("No FP file found in $(folder). Skipping")
        continue
    end

    g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/$(folder)/$(graph_file)", "graphs/$(folder)/$(cfp_file)", "graphs/$(folder)/$(fp_file)")

    gcfp = SafestRoutePair.GraphWithFPandCFP(g, fps, fp_edges, cfps, cfp_edges)

    push!(gcfps, (folder, gcfp))
end



conf_num = 1

# Naive
if run_naive
    println("Naive run started on $(Dates.now())")
    full_folder_name = "$(res_folder_name)/$(conf_num)_naive"

    if !isdir("logs/$(full_folder_name)")
        mkdir("logs/$(full_folder_name)")
    end
    if !isdir("results/$(full_folder_name)")
        mkdir("results/$(full_folder_name)")
    end

    open("results/$(full_folder_name)/params.txt", "w") do io
        println(io, "Naive")
    end

    for (graph_name, gcfp) in gcfps
        println("$(conf_num) $(graph_name) Naive started on $(Dates.now())")
        result = SafestRoutePair.safest_route_pairs_all_naive(gcfp, limit_pairs=limit_pairs,
            use_folds=false)

        results = [result for _ in 1:number_of_runs]

        local_route_df, local_result_df = create_dfs_from_results(results)
        CSV.write(
            "results/$(full_folder_name)/run_routes_$(graph_name).csv",
            local_route_df,
        )
        CSV.write(
            "results/$(full_folder_name)/run_result_$(graph_name).csv",
            local_result_df,
        )
    end

    global conf_num += 1
end

# ACO configurations
for (conf_name, acoS) in configurations_ACO
    # Preparing folders
    println("$(conf_name) started on $(Dates.now())")
    full_folder_name = "$(res_folder_name)/$(conf_num)_$(conf_name)"

    if !isdir("logs/$(full_folder_name)")
        mkdir("logs/$(full_folder_name)")
    end
    if !isdir("results/$(full_folder_name)")
        mkdir("results/$(full_folder_name)")
    end

    open("results/$(full_folder_name)/params.txt", "w") do io
        println(io, "acoS=", acoS)
    end

    for (graph_name, gcfp) in gcfps
        println("$(conf_num): $(graph_name) ACO started on $(Dates.now())")
        results = Folds.map(
            x -> SafestRoutePair.safest_route_pairs_all_aco(
                gcfp;
                acoS,
                logging_file=(x == 1 ? "logs/$(full_folder_name)/$(graph_name)_ACO_$(x).csv" : ""),
                limit_pairs=limit_pairs,
                use_folds=false
            ),
            1:number_of_runs,
        )

        local_route_df, local_result_df = create_dfs_from_results(results)
        CSV.write(
            "results/$(full_folder_name)/run_routes_$(graph_name).csv",
            local_route_df,
        )
        CSV.write(
            "results/$(full_folder_name)/run_result_$(graph_name).csv",
            local_result_df,
        )
    end


    global conf_num += 1
end

# GA configurations
for (conf_name, gaS) in configurations_GA
    # Preparing folders
    println("$(conf_name) started on $(Dates.now())")
    full_folder_name = "$(res_folder_name)/$(conf_num)_$(conf_name)"

    if !isdir("logs/$(full_folder_name)")
        mkdir("logs/$(full_folder_name)")
    end
    if !isdir("results/$(full_folder_name)")
        mkdir("results/$(full_folder_name)")
    end

    open("results/$(full_folder_name)/params.txt", "w") do io
        println(io, "gaS=", gaS)
    end

    for (graph_name, gcfp) in gcfps
        println("$(conf_num): $(graph_name) GA started on $(Dates.now())")
        results = Folds.map(
            x -> SafestRoutePair.safest_route_pairs_all_ga(
                gcfp;
                gaS,
                logging_file=(x == 1 ? "logs/$(full_folder_name)/$(graph_name)_GA_$(x).csv" : ""),
                limit_pairs,
            ),
            1:number_of_runs,
        )

        local_route_df, local_result_df = create_dfs_from_results(results)
        CSV.write(
            "results/$(full_folder_name)/run_routes_$(graph_name).csv",
            local_route_df,
        )
        CSV.write(
            "results/$(full_folder_name)/run_result_$(graph_name).csv",
            local_result_df,
        )
    end


    global conf_num += 1
end


println("Alg_tester_cli ended on $(Dates.now())")

