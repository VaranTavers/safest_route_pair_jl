function preset_seed(popSize, seed)
    [deepcopy(seed) for _ in 1:popSize]
end