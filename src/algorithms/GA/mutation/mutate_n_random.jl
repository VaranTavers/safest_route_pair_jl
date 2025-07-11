function mutate_n_random(partition, n=3)
    vals_to_mutate = randperm(length(partition))
    partition2 = deepcopy(partition)

    for to_mutate in vals_to_mutate[1:n]
        partition2[to_mutate] = rand(0:2)
    end

    partition2
end