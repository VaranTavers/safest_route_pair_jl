using CSV
using DataFrames

function compare_means(df_1, df_2)
    count(df_1[:, :mean] .> df_2[:, :mean])
end

function compare_maxes(df_1, df_2)
    count(df_1[:, :max] .> df_2[:, :max])
end

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

folder_name = "results/2025-04-24_1"

result_prefix = "run_result_"
dfs = []

graphs = Set()

for variation in readdir(folder_name)
    if isdir("$(folder_name)/$(variation)")
        dfs_var = Dict()
        files = readdir("$(folder_name)/$(variation)")

        for file in filter((f) -> occursin(result_prefix, f), files)
            graph_name = file[length(result_prefix)+1:end-4]
            push!(graphs, graph_name)
            dfs_var[graph_name] = CSV.read("$(folder_name)/$(variation)/$(file)", DataFrame)
        end

        push!(dfs, dfs_var)
    end
end

# Simple mean comparison

for graph in graphs
    comp_df = get_comparison_df(dfs, graph, compare_means)
    CSV.write("$(folder_name)/$(graph)_mean_comp.csv", comp_df)
end

# Simple max comparison

for graph in graphs
    comp_df = get_comparison_df(dfs, graph, compare_maxes)
    CSV.write("$(folder_name)/$(graph)_max_comp.csv", comp_df)
end