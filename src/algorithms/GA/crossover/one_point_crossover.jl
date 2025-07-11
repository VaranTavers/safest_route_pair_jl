""" 
Splits the chormosome of each parent at a random location, and the offspring gets the left part from parent 1 and the right one from parent 2.

Example:
```
Parent 1: 0 1 2 2 1 0
Parent 2: 2 1 0 0 1 2

Cut at position 3 (indexed from 1):
Offspring: 0 1 2 0 1 2
```
"""
function one_point_crossover_naive(partition1, partition2)
    point = rand(1:length(partition1))

    res = deepcopy(partition1[1:point])

    append!(res, partition2[point+1:end])

    res
end