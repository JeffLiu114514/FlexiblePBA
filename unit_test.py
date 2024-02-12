from dala import minimal_BER as dala_minimal_BER
from flexible_dala import minimal_BER as flexible_dala_minimal_BER
from flexible_dala import read_from_json, write_to_json, init_model

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


if __name__ == "__main__":
    distributions = init_model()
    eps = 1e-3

    results = {}
    for n in [8, 16, 32]:
        print("----------------------------")
        print(f"Running for {n} levels")
        # print("dala")
        # dala_refined, dala_best_BER = dala_minimal_BER(n, eps, distributions)
        print("flexible_dala")
        flexible_refined, flexible_best_BER = flexible_dala_minimal_BER(n, eps, distributions)
        results[n] = {"flexible_dala": flexible_best_BER} # "dala": dala_best_BER, 
        
    print(results)