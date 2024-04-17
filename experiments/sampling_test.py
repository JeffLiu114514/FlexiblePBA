from dala import minimal_BER as dala_minimal_BER
from flexible_dala import minimal_BER as flexible_dala_minimal_BER
from flexible_dala import read_from_json, write_to_json, init_model
from dala_genmatrix import simulate_error, dump_matrix_sample
import random
import numpy as np
from scipy.stats import ttest_1samp
import matplotlib.pyplot as plt

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
        elif sample_method == "random":
            sampled_distributions[key] = sorted(random.sample(value, sample_size))
    return sampled_distributions

def str_diff(s1, s2):
    assert len(s1) == len(s2)
    diff = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    assert diff <= len(s1)
    return diff

def init_dist():
    gray_coding = \
    {
        4: ["00", "01", "11", "10"],
        8: ["000", "001", "011", "010", "110", "111", "101", "100"],
        16: ["0000", "0001", "0011", "0010", "0110", "0111", "0101", "0100", "1100", "1101", "1111", "1110", "1010", "1011", "1001", "1000"]
    }
    dist_4 = np.zeros((4, 4))
    dist_8 = np.zeros((8, 8))
    dist_16 = np.zeros((16, 16))
    for i in range(0, 4):
        for j in range(0, 4):
            dist_4[i][j] = str_diff(gray_coding[4][i], gray_coding[4][j])
    for i in range(0, 8):
        for j in range(0, 8):
            dist_8[i][j] = str_diff(gray_coding[8][i], gray_coding[8][j])
    for i in range(0, 16):
        for j in range(0, 16):
            dist_16[i][j] = str_diff(gray_coding[16][i], gray_coding[16][j])
    return dist_4, dist_8, dist_16

def report_ber(matrix, n, dist_4, dist_8, dist_16):
    # res = []
    dist = 0
    num_bits = 0
    if n == 4:
        dist = dist_4
        num_bits = 2
    elif n == 8:
        dist = dist_8
        num_bits = 3
    elif n == 16:
        dist = dist_16
        num_bits = 4
    else:
        assert False
            
    ber_matrix = np.multiply(matrix, dist) / n
    # pprint.pprint(ber_matrix)
    ber_avg = np.sum(ber_matrix) / num_bits
    # res.append(ber_avg)
    return ber_avg

if __name__ == "__main__":
    distributions = init_model()
    eps = 1e-3
    outfile = "./ember_capacity/"
    results = {}
    sample_method="random" #"uniform" or "random"
    num_trails = 50
    dist_4, dist_8, dist_16 = init_dist()
    PRINT_LOG = False
    
    # print("----------------------------")
    # print(f"Running for 4 levels")
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
        
        ber_reduction_list = {}
        stats_significance = {}
        
        for sample_size in range(15, 100):
            dala_ber_list = []
            fdala_ber_list = []
            for i in range(num_trails):
                try:
                    sample_distribution = sample_distributions(distributions, sample_size, sample_method=sample_method)
                    if PRINT_LOG: print(f"Sample size: {sample_size}%")
                    
                    if PRINT_LOG: print("dala")
                    dala_refined, dala_best_BER = dala_minimal_BER(n, eps, sample_distribution)
                    dala_P = simulate_error(dala_refined, distributions)
                    # dump_matrix_sample(P, outfile + "dala", sample_size, suffix=sample_method)
                    dala_ber = report_ber(dala_P, n, dist_4, dist_8, dist_16)
                    
                    if PRINT_LOG: print("flexible_dala")
                    flexible_refined, flexible_best_BER, _ = flexible_dala_minimal_BER(n, eps, sample_distribution)
                    fdala_P = simulate_error(flexible_refined, distributions)
                    # dump_matrix_sample(P, outfile + "flexible", sample_size, suffix=sample_method)
                    fdala_ber = report_ber(fdala_P, n, dist_4, dist_8, dist_16)
                    
                    dala_ber_list.append(dala_ber)
                    fdala_ber_list.append(fdala_ber)
                    
                    # fkey = f"{n}_{sample_size}"
                    # results[(n, sample_size)] = {"dala": (dala_refined, dala_best_BER), "flexible_dala": (flexible_refined, flexible_best_BER)} # (n, sample_size)
                except UnboundLocalError as e:
                    if PRINT_LOG: print(f"Sample size: {sample_size}%")
                    if PRINT_LOG: print(f"Error: {e}")
                    continue
            
            if len(dala_ber_list) == 0 or len(fdala_ber_list) == 0:
                continue
            assert len(dala_ber_list) == len(fdala_ber_list)
            
            
            ber_reduction = [(dala_ber_list[j] - fdala_ber_list[j]) for j in range(len(dala_ber_list))]# / dala_ber_list[j]
            t_stat, p_value = ttest_1samp(ber_reduction, 0, alternative='greater')
            stats_significance[sample_size] = p_value
            # print(dala_ber_list)
            # print(fdala_ber_list)
            # print(ber_reduction)
            ber_reduction_list[sample_size] = [np.mean(ber_reduction), np.std(ber_reduction)]
        
        results[n] = ber_reduction_list
        
        ber_reduction_values = []
        for key, value in ber_reduction_list.items():
            if PRINT_LOG: print(key, value)
            if PRINT_LOG: print()
            ber_reduction_values.append(value[0])
        
        t_stat, p_value = ttest_1samp(ber_reduction_values, 0, alternative='greater') # Change the condition to 'greater'
        significance = "significant" if p_value < 0.05 else "not significant"
        print(f"Average BER reduction for {n} levels sampling with {num_trails} trails has p value of {p_value} and is {significance}")
        for key, value in stats_significance.items():
            if value < 0.05:
                print(f"Sample size {key} has p value of {value} and is significant")
            else:
                print(f"Sample size {key} has p value of {value} and is not significant")
        
        
        # plt.errorbar(ber_reduction_list.keys(), [v[0] for v in ber_reduction_list.values()], yerr=[v[1] for v in ber_reduction_list.values()], marker='o')
        # # plt.plot(ber_reduction_list.keys(), ber_reduction_list.values(), marker='o')
        # plt.xlabel('sample sizes')
        # plt.ylabel('BER reductions')
        # plt.title(f'{n}-level BER reductions')
        # plt.show()
        
        # plt.plot()
            
    
    
    # write_to_json(results, "sampling_test_results.json")
    

