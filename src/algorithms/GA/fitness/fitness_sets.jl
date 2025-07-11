# Alternative fitness function: sum of FP-s from both sets should be maximal
function calc_fitness_sets(solution, g_indep, runS::GaRunSettings)  #_simple

    paths = partition_to_paths(g_indep, solution, runS.fp_edges, runS.source, runS.target)


    if fst(paths) == [] && snd(paths) == []
        return -20
    elseif fst(paths) == []
        return -10
    elseif snd(paths) == []
        return -10
    end

    solution2 = deepcopy(solution)

    if length(runS.fp_depend) === length(runS.fp_edges)
        edge_comps = runS.fp_depend

        for (i, (sol, comp)) in enumerate(zip(solution, edge_comps))
            l = length(comp)
            if l > 1 && sol == 0
                vals = [solution[j] for j in comp]
                if all(vals .== 1)
                    solution2[i] = 1
                elseif all(vals .== 2)
                    solution2[i] = 2
                end
            end
        end
    end

    sum(runS.fps[solution2.!=0])
end
