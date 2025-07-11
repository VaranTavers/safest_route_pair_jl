function mutate_one_random(partition)
    to_mutate = rand(1:length(partition))
    partition2 = deepcopy(partition)
    partition2[to_mutate] = rand(0:2)

    partition2
end