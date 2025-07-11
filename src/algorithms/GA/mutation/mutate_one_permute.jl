# Basic mutation that permutates a value in the 0->1->2->0 order on one of the genes
function mutate_one_permute(partition)
    to_mutate = rand(1:length(partition))
    partition2 = deepcopy(partition)
    partition2[to_mutate] = (partition[to_mutate] + 1) % 3

    partition2
end