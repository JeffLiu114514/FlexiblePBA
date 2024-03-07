from dala import minimal_BER as dala_minimal_BER
from flexible_dala import minimal_BER as flexible_dala_minimal_BER
from flexible_dala import read_from_json, write_to_json, init_model
from dala_genmatrix import simulate_error, dump_matrix_sample
import random

def dump_to_json(level_alloc, algo_name):
    if len(level_alloc) == 16:
        bits_per_cell = 4
    elif len(level_alloc) == 8:
        bits_per_cell = 3
    elif len(level_alloc) == 4:
        bits_per_cell = 2
    elif len(level_alloc) == 2:
        bits_per_cell = 1
    elif len(level_alloc) == 32:
        bits_per_cell = 5
    bpc = read_from_json(f"settings/{bits_per_cell}bpc.json")
    for i in range(0, len(level_alloc)):
        # [Rlow, Rhigh, tmin, tmax]
        bpc['level_settings'][i]["adc_upper_read_ref_lvl"] = level_alloc[i][1]
        bpc['level_settings'][i]["adc_lower_write_ref_lvl"] = level_alloc[i][2]
        bpc['level_settings'][i]["adc_upper_write_ref_lvl"] = level_alloc[i][3]
    write_to_json(bpc, f"unit_tests/{bits_per_cell}bpc_{algo_name}.json")
    

def sample_distributions(distributions, percentage, sample_method="uniform"):
    sampled_distributions = {}
    for key, value in distributions.items():
        sample_size = int(len(value) * (percentage / 100))
        if sample_method == "uniform":
            sampled_distributions[key] = value[::len(value)//sample_size]
        else:
            sampled_distributions[key] = random.sample(value, sample_size)
    return sampled_distributions


if __name__ == "__main__":
    distributions = init_model()
    eps = 1e-3
    outfile = "./ember_capacity/"
    results = {}
    
    print("----------------------------")
    print(f"Running for 4 levels")
    # for sample_size in [25, 50, 75]:
    #     n = 4
    #     sample_distribution = sample_distributions(distributions, sample_size)
        
    #     print("dala")
    #     dala_refined, dala_best_BER = dala_minimal_BER(n, eps, sample_distribution, 0, 1, True)
    #     P = simulate_error(dala_refined, distributions)
    #     dump_matrix_sample(P, outfile + "dala", sample_size)
    #     print("flexible_dala")
    #     flexible_refined, flexible_best_BER = flexible_dala_minimal_BER(n, eps, sample_distribution, 0, 1, True)
    #     P = simulate_error(flexible_refined, distributions)
    #     dump_matrix_sample(P, outfile + "flexible", sample_size)
        
    #     results[(n, sample_size)] = {"dala": (dala_refined, dala_best_BER), "flexible_dala": (flexible_refined, flexible_best_BER)}
    
    for n in [8, 16]:
        print("----------------------------")
        print(f"Running for {n} levels")
        
        for sample_size in range(15, 100):
            try:
                sample_distribution = sample_distributions(distributions, sample_size)
                print(f"Sample size: {sample_size}%")
                print("dala")
                dala_refined, dala_best_BER = dala_minimal_BER(n, eps, sample_distribution)
                P = simulate_error(dala_refined, distributions)
                dump_matrix_sample(P, outfile + "dala", sample_size)
                print("flexible_dala")
                flexible_refined, flexible_best_BER = flexible_dala_minimal_BER(n, eps, sample_distribution)
                P = simulate_error(flexible_refined, distributions)
                dump_matrix_sample(P, outfile + "flexible", sample_size)
                # fkey = f"{n}_{sample_size}"
                results[(n, sample_size)] = {"dala": (dala_refined, dala_best_BER), "flexible_dala": (flexible_refined, flexible_best_BER)} # (n, sample_size)
            except UnboundLocalError as e:
                print(f"Sample size: {sample_size}%")
                print(f"Error: {e}")
                continue
    
    # for key, value in results.items():
    #     print(key, value)
    #     print()
    # write_to_json(results, "sampling_test_results.json")
    

