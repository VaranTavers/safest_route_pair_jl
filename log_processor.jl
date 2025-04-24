using CSV
using DataFrames
using Plots

function vals_to_vec(type, str)
    if startswith(str, "Real")
        str = str[5:end]
    end
    @assert str[1] == '[' && str[end] == ']'

    valsStr = split(replace(str, '[' => "", ']' => ""), ',')

    [parse(type, val) for val in valsStr]
end

folder_name = "logs/ACO_DEFAULT_GA_SIMPLER_GENES_2025-04-23_3/"

folder = readdir(folder_name)

for file in folder
    if occursin("GA", file)
        df = CSV.read("$(folder_name)/$(file)", DataFrames.DataFrame, header=false)

        vals = [vals_to_vec(Float32, row["Column3"]) for row in eachrow(df)]


        A = [x[1] for x in vals]
        #@show A
        B = [maximum(A[1:i]) for i in eachindex(A)]
        #@show B

        @show file, count(B .< B[end])
    end
end

