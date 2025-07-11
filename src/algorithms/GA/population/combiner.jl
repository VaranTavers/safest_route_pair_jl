function combine_seed(f1, f2, cSize, popSize, rate)
    n1 = Int(floor(popSize * rate))
    r1 = f1(n1, cSize)
    r2 = f2(popSize - n1, cSize)

    append!(r1, r2)

    r1
end