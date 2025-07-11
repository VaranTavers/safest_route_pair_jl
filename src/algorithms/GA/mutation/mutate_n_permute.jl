# Basic mutation that permutates n values in the 0->1->2->0 order on one of the genes
function mutate_n_permute(partition, n=3)
    vals_to_mutate = randperm(length(partition))
    partition2 = deepcopy(partition)

    for to_mutate in vals_to_mutate[1:n]
        partition2[to_mutate] = (partition[to_mutate] + 1) % 3
    end

    partition2
end