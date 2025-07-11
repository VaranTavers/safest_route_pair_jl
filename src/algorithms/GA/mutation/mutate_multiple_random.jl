function mutate_multiple_random(partition, gene_mut_prob=0.01)
    [rand() < gene_mut_prob ? rand(0:2) : x for x in partition]
end