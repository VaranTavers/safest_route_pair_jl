# Basic crossover that chooses each cromosome randomly from one of the parents
function npoint_crossover_naive(partition1, partition2)
    [rand((1, 2)) == 1 ? partition1[i] : partition2[i] for i in eachindex(partition1)]
end