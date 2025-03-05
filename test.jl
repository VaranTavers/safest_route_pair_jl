using Pkg

Pkg.activate(".")

using SafestRoutePair
using Graphs
using DataFrames
using CSV
using Statistics
using Dates
using Folds


#=g, probs, edges = AntSafestRoute.read_graph_and_failure(
    "../interroute_j.gml",
    "interroute_v2",
    "../whole_graph_complete_IX_upper_big_grid.xml",
)


AntSafestRoute.safestRoutesAll(
    "../interroute_j.gml",
    "interroute_v2",
    "../whole_graph_complete_IX_upper_big_grid.xml",
    "GML+XML",
)=#

folders = readdir("./graphs")


println("Alg_tester_cli started on $(Dates.now())")

fst((x, _)) = x
snd((_, y)) = y

function toString(x)
    "$(x)"
end

numberOfRuns = 1

param_tuning_α = [1]
param_tuning_β = [1.5]
param_tuning_nr_ants = [25]
param_tuning_ρ = [0.3]
param_tuning_ϵ = [0.1]
param_tuning_nr_gen = [1]
param_tuning_starting_pheromone = [1]

configurations = [
    # conf_name,                        n_p, mut, cro, elit, crossoverAlg, mutationAlg,       meme,  log,  iter 
    # Baselines
    ("ACO_DEFAULT", ACOSettings(α, β, nr_ants, ρ, ϵ, nr_gen, s_pheromone)) for
    α in param_tuning_α, β in param_tuning_β, nr_ants in param_tuning_nr_ants,
    ρ in param_tuning_ρ, ϵ in param_tuning_ϵ, nr_gen in param_tuning_nr_gen,
    s_pheromone in param_tuning_starting_pheromone
]

if !isdir("./logs")
    mkdir("./logs")
end
if !isdir("./results")
    mkdir("./results")
end

for (ii, (conf_name, acoS)) in enumerate(configurations)
    println("$(conf_name) started on $(Dates.now())")

    date_of_start = Dates.today()

    res_folder_name = "$(conf_name)_$(ii)_$(date_of_start)"
    if !isdir("logs/$(res_folder_name)")
        mkdir("logs/$(res_folder_name)")
    end
    if !isdir("results/$(res_folder_name)")
        mkdir("results/$(res_folder_name)")
    end


    open("results/$(res_folder_name)/params.txt", "a") do io
        println(io, "acoS=", acoS)
    end

    for (i, folder) in enumerate(folders)
        files = readdir("graphs/$(folder)")
        println("$(i)/$(length(folders)) $(folder) started on $(Dates.now())")

        graph_files = collect(filter(x -> x[(end-3):end] == ".gml", files))
        if length(graph_files) < 1
            println("No GML file found in $(folder). Skipping")
            continue
        end

        graph_file = graph_files[1]
        cfp_file = "cfp.xml"
        fp_file = "fp.xml"


        local_route_df = DataFrame()
        local_result_df = DataFrame()

        g, cfps, cfp_edges, fps, fp_edges = read_graph_and_failure("graphs/$(folder)/$(graph_file)", "graphs/$(folder)/$(cfp_file)", "graphs/$(folder)/$(fp_file)")

        gcfp = SafestRoutePair.GraphWithFPandCFP(g, fps, fp_edges, cfps, cfp_edges)

        results = Folds.map(
            x -> SafestRoutePair.safest_route_pairs_all_aco(
                gcfp;
                acoS,
                logging_file=(x == 1 ? "logs/$(res_folder_name)/$(folder)_$(x).csv" : "")
            ),
            1:numberOfRuns,
        )

        local_route_df[!, :nodes] = fst.(results[1])
        local_result_df[!, :nodes] = fst.(results[1])

        for (i, result) in enumerate(results)
            local_route_df[!, :x] = toString.(fst.(snd.(result)))
            rename!(local_route_df, :x => "aco$(i)")
            local_result_df[!, :x] = snd.(snd.(result))
            rename!(local_result_df, :x => "aco$(i)")
        end

        tmpMat = Matrix(local_result_df[:, 2:(numberOfRuns+1)])
        means = mean(tmpMat, dims=2)
        stds = std(tmpMat, dims=2)
        mins = minimum(tmpMat, dims=2)
        maxs = maximum(tmpMat, dims=2)


        local_result_df[!, :mean_aco] = means[:]
        local_result_df[!, :std_aco] = stds[:]
        local_result_df[!, :min_aco] = mins[:]
        local_result_df[!, :max_aco] = maxs[:]

        results2 = Folds.map(
            x -> SafestRoutePair.safest_route_pairs_all_ga(
                gcfp;
                logging_file=(x == 1 ? "logs/$(res_folder_name)/$(folder)_$(x).csv" : "")
            ),
            1:numberOfRuns,
        )


        for (i, result) in enumerate(results2)
            local_route_df[!, :x] = toString.(fst.(snd.(result)))
            rename!(local_route_df, :x => "ga$(i)")
            local_result_df[!, :x] = snd.(snd.(result))
            rename!(local_result_df, :x => "ga$(i)")
        end

        tmpMat = Matrix(local_result_df[:, 2:(numberOfRuns+1)])
        means = mean(tmpMat, dims=2)
        stds = std(tmpMat, dims=2)
        mins = minimum(tmpMat, dims=2)
        maxs = maximum(tmpMat, dims=2)


        local_result_df[!, :mean_ga] = means[:]
        local_result_df[!, :std_ga] = stds[:]
        local_result_df[!, :min_ga] = mins[:]
        local_result_df[!, :max_ga] = maxs[:]



        CSV.write(
            "results/$(res_folder_name)/run_routes_$(folder).csv",
            local_route_df,
        )
        CSV.write(
            "results/$(res_folder_name)/run_result_$(folder).csv",
            local_result_df,
        )
    end
end


println("Alg_tester_cli ended on $(Dates.now())")

