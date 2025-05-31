# General configuration
number_of_runs = 20
limit_pairs = 0
run_naive = true
run_suurballe = true

# ACO configuration
aco_param_tuning_nr_gen = [5, 100]
aco_param_tuning_α = [1]
aco_param_tuning_β = [1.5]
aco_param_tuning_nr_ants = [25]
aco_param_tuning_ρ = [0.3]
aco_param_tuning_ϵ = [0.1]
aco_param_tuning_starting_pheromone = [1]

# GA configuration
ga_param_tuning_nr_gen = [100]
ga_param_tuning_n_p = [50]
ga_param_tuning_mut_rate = [0.5]
ga_param_tuning_cro_rate = [0.9]
ga_param_tuning_elit = [0.5]
ga_param_tuning_crossover = [npoint_crossover_naive, one_point_crossover_naive]
ga_param_tuning_mutation = [mutate_permute, mutate_random]
ga_param_tuning_fitness = [calc_fitness_paths, calc_fitness_sets] # 
ga_param_prob_limit = [0]
ga_param_edge_limit = [0]
ga_param_seed_naive = [false]