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

def decide_end_level(point, level_alloc):
    for i in range(len(level_alloc)):
        rmin, rmax, wmin, wmax = level_alloc[i]    
        if point >= rmin and point < rmax:
            return i
    if point == 64: # last level
        return len(level_alloc)-1
    print(level_alloc)
    assert False, point
    
def simulate_error(level_alloc, distributions):
    num_levels = len(level_alloc)
    P = np.zeros((num_levels, num_levels))
    num_points = 0
    for i in range(len(level_alloc)):
        rmin, rmax, wmin, wmax = level_alloc[i]
        d1s = distributions[(wmin, wmax)]
        for point in d1s:
            end_level = decide_end_level(point, level_alloc)
            P[end_level][i] += 1
        for j in range(0, len(level_alloc)):
            P[j][i] = P[j][i] / len(d1s)
        num_points += len(d1s)
    return P

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