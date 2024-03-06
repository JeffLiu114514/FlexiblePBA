import itertools
import numpy as np
import time

gray_coding = \
{
    4: ["00", "01", "11", "10"],
    8: ["000", "001", "011", "010", "110", "111", "101", "100"],
    16: ["0000", "0001", "0011", "0010", "0110", "0111", "0101", "0100", "1100", "1101", "1111", "1110", "1010", "1011", "1001", "1000"]
}

def str_diff(s1, s2):
    assert len(s1) == len(s2)
    diff = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    assert diff <= len(s1)
    return diff

def get_matrix_from_file(filename):
    with open(filename, "r") as fin:
        lines = fin.readlines()
        n = len(lines)
        matrix = np.zeros((n, n))
        for i in range(n):
            line = list((map(float, lines[i].split(","))))
            assert len(line) == n
            for j in range(n):
                matrix[i][j] = line[j]
    return matrix

def get_best(file_prefix, num):
    fname = file_prefix + str(num)
    matrix = get_matrix_from_file(fname)
    best_ber = 9999999999
    best_coding = []
    cnt = 0
    time_start = time.time()
    for p in itertools.permutations(gray_coding[num]):
        dist = np.zeros((num, num))
        for i in range(0, num):
            for j in range(0, num):
                dist[i][j] = str_diff(p[i], p[j])
        ber_matrix = np.multiply(matrix, dist) / num
        ber_avg = np.sum(ber_matrix) / np.log2(num)
        if ber_avg < best_ber:
            best_ber = ber_avg
            best_coding = p
            
        cnt += 1
        if cnt % 100000 == 0:
            print(cnt)
            time_cut = time.time()
            print(time_cut - time_start)
            print("estimated time: ", (time_cut - time_start) * 2*10**13 / cnt)
 
    print(best_ber, best_coding)
    
    return best_ber, best_coding

if __name__ == "__main__":
    # get_best("dala", 4)
    # get_best("dala", 8)
    get_best("dala", 16)
    
    # get_best("flexible", 4)
    # get_best("flexible", 8)
    # get_best("flexible", 16)


