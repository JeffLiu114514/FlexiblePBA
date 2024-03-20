from dala import minimal_BER as dala_minimal_BER
from flexible_dala import minimal_BER as flexible_dala_minimal_BER
from flexible_dala import read_from_json, write_to_json, init_model
from dala_genmatrix import simulate_error
import numpy as np


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

def get_ber_for_allocs(dala_alloc, distributions, n, dist_4, dist_8, dist_16):
    # dist_4, dist_8, dist_16 = init_dist()
    P = simulate_error(dala_alloc, distributions)
    ber = report_ber(P, n, dist_4, dist_8, dist_16)
    return ber

if __name__ == "__main__":
    distributions = init_model()
    eps = 1e-5

    results = {}
    level_list = [8, 16]
    for n in level_list:
        # get level alloc
        print("----------------------------")
        print(f"Running for {n} levels")
        print("dala")
        dala_refined, dala_best_BER = dala_minimal_BER(n, eps, distributions)
        print("flexible_dala")
        flexible_refined, flexible_best_BER = flexible_dala_minimal_BER(n, eps, distributions)
        
        # simulate error
        dala_P = simulate_error(dala_refined, distributions)
        fdala_P = simulate_error(flexible_refined, distributions)
        
        # report ber
        dist_4, dist_8, dist_16 = init_dist()
        dala_ber = report_ber(dala_P, n, dist_4, dist_8, dist_16)
        fdala_ber = report_ber(fdala_P, n, dist_4, dist_8, dist_16)
        
        results[n] = {"fdala_ber": fdala_ber, "dala_ber": dala_ber}
        
    print(results)