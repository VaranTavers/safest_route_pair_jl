function random_seed(cSize, popSize)
    [[rand([1, 2]) for _ in 1:cSize] for _ in 1:popSize]
end