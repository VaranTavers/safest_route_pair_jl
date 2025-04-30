using CSV
using DataFrames
using Plots
using HypothesisTests
using Statistics
using StatsPlots

# Simple comparisons

function compare_means(df_1, df_2)
    count(df_1[:, :mean] .> df_2[:, :mean])
end

function compare_maxes(df_1, df_2)
    count(df_1[:, :max] .> df_2[:, :max])
end

# Chess ranking functions

function g(Φ)
    1 / sqrt(1 + 3 * Φ^2 / π^2)
end

function E(μ, μ_j, Φ_j)
    1 / (1 + exp(-g(Φ_j) * (μ - μ_j)))
end

function compare_random(df_1, df_2, row)
    _r1, c1 = size(df_1)
    _r2, c2 = size(df_2)
    o1 = rand(2:(c1-4))
    o2 = rand(2:(c2-4))

    #@show row, df_1[row, 1], df_2[row, 1], df_1[row, o1], df_2[row, o2]
    if df_1[row, o1] == df_2[row, o2]
        return 0.5
    end
    if df_1[row, o1] > df_2[row, o2]
        return 1
    end

    0
end


function get_pairings(dfs, graph, num_challengers, num_challenges)
    @assert num_challenges <= num_challengers - 1
    competitors = [[] for _ in 1:num_challengers]
    results = [[] for _ in 1:num_challengers]

    for i in 1:num_challengers
        while length(competitors[i]) < num_challenges
            new = rand(1:num_challengers)
            while new == i || new in competitors[i]
                new = rand(1:num_challengers)
            end

            #@show i, new
            possible_rows, _ = size(dfs[1][graph])
            result = compare_random(dfs[i][graph], dfs[new][graph], rand(1:possible_rows))
            push!(competitors[i], new)
            push!(competitors[new], i)
            push!(results[i], result)
            push!(results[new], 1 - result)
        end
    end

    competitors, results
end

function f(x, Δ, Φ, v, a, τ)
    exp(x) * (Δ^2 - Φ^2 - v - exp(x)) / 2(Φ^2 + v + exp(x))^2 - (x - a) / τ^2
end

function calculate_new_sigma(Δ, Φ, v, τ, σ)
    ϵ = 0.000001
    a = log(σ^2)
    A = a
    B = 0

    local_f(x) = f(x, Δ, Φ, v, a, τ)
    if Δ^2 > Φ^2 + v
        B = log(Δ^2 - Φ^2 - v)
    else
        k = 1
        while local_f(a - k * τ) < 0
            k += 1
        end
        B = a - k * τ
    end

    fA = local_f(A)
    fB = local_f(B)
    while abs(B - A) > ϵ
        C = A + (A - B) * fA / (fB - fA)
        fC = local_f(C)
        if fB * fC <= 0
            A = B
            fA = fB
        else
            fA /= 2
        end
        B = C
        fB = fC
    end

    exp(A / 2)
end

# Based on: https://www.sciencedirect.com/science/article/pii/S002002551400276X
# and:      https://www.glicko.net/glicko/glicko2.pdf
function chess_ranking(dfs, num_tournamets, num_challenges)
    num_competitors = length(dfs)

    τ = 0.9 # Constraints the volatility
    ratings = 1500 .* ones(num_competitors)
    rating_deviations = 350 .* ones(num_competitors)
    rating_volatilities = 0.06 .* ones(num_competitors) # σ

    for i in 1:num_tournamets
        μ = (ratings .- 1500) ./ 173.7178
        Φ = rating_deviations ./ 173.7178

        graph_to_test = rand(keys(dfs[1]))
        #@show graph_to_test
        pairings, s_s = get_pairings(dfs, graph_to_test, num_competitors, num_challenges)



        #@show pairings, s_s
        v = [sum([g(Φ[j])^2 * E(μ[i], μ[j], Φ[j]) * (1 - E(μ[i], μ[j], Φ[j])) for j in pairings[i]])^-1 for i in 1:num_competitors]


        Δ = [v[i] * sum([g(Φ[j]) * (s_s[i][k] - E(μ[i], μ[j], Φ[j])) for (k, j) in enumerate(pairings[i])]) for i in 1:num_competitors]


        rating_volatilities = [calculate_new_sigma(Δ[i], Φ[i], v[i], τ, rating_volatilities[i]) for i in 1:num_competitors]
        Φ_pre = [sqrt(Φ_i^2 + σ^2) for (Φ_i, σ) in zip(Φ, rating_volatilities)]

        Φ = [1 / sqrt(1 / Φ_i^2 + 1 / v[i]) for (i, Φ_i) in enumerate(Φ_pre)]
        μ = μ + [Φ[i]^2 * sum([g(Φ[j]) * (s_s[i][k] - E(μ[i], μ[j], Φ[j])) for (k, j) in enumerate(pairings[i])]) for i in 1:num_competitors]

        ratings = 173.7178 .* μ .+ 1500
        rating_deviations = 173.7178 .* Φ
    end

    ratings, rating_deviations
end


# General funtions

function get_comparison_df(dfs, graph, comp_f)
    l = length(dfs)
    mat = zeros(l, l)

    for i in 1:l
        for j in 1:l
            if i != j
                mat[i, j] = comp_f(dfs[i][graph], dfs[j][graph])
            end
        end
    end

    comp_df = DataFrame()
    comp_df[!, :nodes] = collect(1:l)

    for (i, val) in enumerate(eachcol(mat))
        comp_df[!, :x] = val
        rename!(comp_df, :x => "$(i)")
    end

    comp_df
end

### MAIN PROGRAM

run_mean_comp = true
run_max_comp = true
run_chess = true
run_wilcoxon = true
number_of_heatmaps_per_graph = 5

folder_name = "results/2025-04-27_2"

result_prefix = "run_result_"
dfs = []

graphs = Set()

variations = collect(filter(x -> (isdir("$(folder_name)/$(x)") && x != "images"), readdir(folder_name)))
variation_nums = [parse(Int, split(x, "_")[1]) for x in variations]
variations = variations[sortperm(variation_nums)]
@show variations


for variation in variations
    dfs_var = Dict()
    files = readdir("$(folder_name)/$(variation)")

    for file in filter((f) -> occursin(result_prefix, f), files)
        graph_name = file[length(result_prefix)+1:end-4]
        push!(graphs, graph_name)
        dfs_var[graph_name] = CSV.read("$(folder_name)/$(variation)/$(file)", DataFrame)
    end

    push!(dfs, dfs_var)
end

# Simple mean comparison
if run_mean_comp
    for graph in graphs
        comp_df = get_comparison_df(dfs, graph, compare_means)
        CSV.write("$(folder_name)/$(graph)_mean_comp.csv", comp_df)
    end
end

# Simple max comparison

if run_max_comp
    for graph in graphs
        comp_df = get_comparison_df(dfs, graph, compare_maxes)
        CSV.write("$(folder_name)/$(graph)_max_comp.csv", comp_df)
    end
end

if run_chess
    chess_result, confidence = chess_ranking(dfs, 10000, length(dfs) - 1)

    @show chess_result, confidence

    sorted_indices = sortperm(chess_result, rev=true)
    sorted_result = chess_result[sorted_indices]
    sorted_confidence = confidence[sorted_indices]

    chess_p = scatter(sorted_result, [10 * i for i in length(chess_result):-1:0], yticks=([10 * i for i in length(chess_result):-1:0], sorted_indices), xerror=sorted_confidence, title="Chess ranking", label=false)

    savefig(chess_p, "$(folder_name)/chess_res.pdf")
end

if run_wilcoxon
    if !isdir("$(folder_name)/images")
        mkdir("$(folder_name)/images")
    end
    number_of_variations = length(dfs)

    for graph in graphs
        point_sum = zeros(number_of_variations, number_of_variations)
        point_sum_d = zeros(number_of_variations, number_of_variations)

        l_df, _ = size(dfs[1][graph])

        if !isdir("$(folder_name)/images/$(graph)")
            mkdir("$(folder_name)/images/$(graph)")
        end

        for node_pair in 1:l_df
            conf_mat = [
                pvalue(MannWhitneyUTest(Vector(df_i[graph][node_pair, 2:end-4]), Vector(df_j[graph][node_pair, 2:end-4]))) < 0.05 ?
                (mean(df_i[graph][node_pair, :mean]) > mean(df_j[graph][node_pair, :mean]) ? 1 : -1) : 0 for
                df_i in dfs, df_j in dfs
            ]
            conf_mat_d = [
                pvalue(MannWhitneyUTest(Vector(df_i[graph][node_pair, 2:end-4]), Vector(df_j[graph][node_pair, 2:end-4]))) < 0.05 ? 1 : 0 for
                df_i in dfs, df_j in dfs
            ]

            df = DataFrame(x=Int[], y=Float64[])

            for (i, df_i) in enumerate(dfs)
                for v in df_i[graph][node_pair, 2:end-4]
                    push!(df, (i, v))
                end
            end


            point_sum += conf_mat
            point_sum_d += conf_mat_d

            if node_pair <= number_of_heatmaps_per_graph
                savefig(
                    Plots.heatmap(
                        conf_mat_d,
                        yflip=true,
                        aspect_ratio=0.5,
                        color=[:white, :black],
                        showaxis=:xy,
                        xticks=1:number_of_variations,
                        yticks=1:number_of_variations,
                        title="$(graph)_$(dfs[1][graph][node_pair, 1])",
                    ),
                    "$(folder_name)/images/$(graph)/heat_$(graph)_$(dfs[1][graph][node_pair, 1]).pdf",
                )
                savefig(
                    Plots.heatmap(
                        conf_mat,
                        yflip=true,
                        aspect_ratio=0.5,
                        showaxis=:xy,
                        xticks=1:number_of_variations,
                        yticks=1:number_of_variations,
                        title="$(graph)_$(dfs[1][graph][node_pair, 1])",
                    ),
                    "$(folder_name)/images/$(graph)/heat2_$(graph)_$(dfs[1][graph][node_pair, 1]).pdf",
                )

                savefig(
                    boxplot(
                        df[:, :x],
                        df[:, :y],
                        line=(2, :black),
                        fill=(0.3, :orange),
                        legend=false,
                        xticks=1:number_of_variations,
                        ylabel="fitness value",
                        title="$(graph)_$(dfs[1][graph][node_pair, 1])",
                    ),
                    "$(folder_name)/images/$(graph)/box_$(graph)_$(dfs[1][graph][node_pair, 1]).pdf",
                )



                savefig(
                    scatter(
                        df[:, :x],
                        df[:, :y],
                        legend=false,
                        xticks=1:number_of_variations,
                        title="$(graph)_$(dfs[1][graph][node_pair, 1])",
                    ),
                    "$(folder_name)/images/$(graph)/scatter_$(graph)_$(dfs[1][graph][node_pair, 1]).pdf",
                )

                df = DataFrame(conf_mat, ["$(i)" for i in 1:number_of_variations])
                df[!, "Config"] = ["$(i)" for i in 1:number_of_variations]
                CSV.write("$(folder_name)/images/$(graph)/$(graph)_$(dfs[1][graph][node_pair, 1]).csv", df)
            end
        end
        point_sum
        savefig(
            heatmap(
                point_sum,
                yflip=true,
                aspect_ratio=0.5,
                showaxis=:xy,
                xticks=1:number_of_variations,
                yticks=1:number_of_variations,
            ),
            "$(folder_name)/images/sum_$(graph).pdf",
        )

        point_sum_d ./= l_df
        savefig(
            heatmap(
                point_sum_d,
                yflip=true,
                aspect_ratio=0.5,
                color=[:white, :black],
                showaxis=:xy,
                xticks=1:number_of_variations,
                yticks=1:number_of_variations,
            ),
            "$(folder_name)/images/percent_$(graph).pdf",
        )

        sum_df = DataFrame(point_sum, ["$(i)" for i in 1:number_of_variations])
        sum_df[!, "Config"] = ["$(i)" for i in 1:number_of_variations]

        sum_of_rows = sum(point_sum, dims=2)
        best_version = argmax(sum_of_rows)[1]
        best_value = maximum(sum_of_rows)
        @show sum_of_rows
        sum_df[!, "Sum"] = sum_of_rows[:]

        CSV.write("$(folder_name)/images/sum_$(graph).csv", sum_df)


    end

end