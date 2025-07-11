function random_seed_only_sets(cSize, popSize)
    [[rand([1, 2]) for _ in 1:cSize] for _ in 1:popSize]
end