import numpy as np
import pprint
import matplotlib.pyplot as plt

directory = "./ember_capacity/"

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

# def write_matrix(matrix, filename):
#     with open(filename, "w") as fout:
#         to_write = []
#         for i in range(len(matrix)):
#             to_write.append(",".join(map(str, matrix[i])) + "\n")
#         fout.writelines(to_write)

def compute_average(matrix):
    error_rate = 0
    for i in range(len(matrix)):
        error_rate += 1 - matrix[i][i]
    error_rate = error_rate / len(matrix)
    return error_rate
            
def report_results(filename_prefix, hint):
    res = {}
    for i in range(4, 9):
        fname = filename_prefix + str(i)
        matrix = get_matrix_from_file(fname)
        avg = compute_average(matrix)
        res[i] = avg
    print(hint + " = \\")
    pprint.pprint(res)
    return res

def report_drift_reduction(queries):
    for item in queries:
        our, pba, hint = item
        assert len(our) == len(pba)
        reduce_list = []
        for i in range(4, 17):
            reduce_list.append((pba[i] - our[i]) / pba[i])
        reduce_avg = sum(reduce_list) / len(reduce_list)
        print(f"{hint} Drift Reduction", reduce_list)
        print(f"{hint} Average Drift Reduction", reduce_avg)

gray_coding = \
{
    4: ["00", "01", "11", "10"],
    8: ["000", "001", "011", "010", "110", "111", "101", "100"],
    16: ["0000", "0001", "0011", "0010", "0110", "0111", "0101", "0100", "1100", "1101", "1111", "1110", "1010", "1011", "1001", "1000"]
}
dist_4 = np.zeros((4, 4))
dist_8 = np.zeros((8, 8))
dist_16 = np.zeros((16, 16))

def str_diff(s1, s2):
    assert len(s1) == len(s2)
    diff = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    assert diff <= len(s1)
    return diff

def init_dist():
    for i in range(0, 4):
        for j in range(0, 4):
            dist_4[i][j] = str_diff(gray_coding[4][i], gray_coding[4][j])
    for i in range(0, 8):
        for j in range(0, 8):
            dist_8[i][j] = str_diff(gray_coding[8][i], gray_coding[8][j])
    for i in range(0, 16):
        for j in range(0, 16):
            dist_16[i][j] = str_diff(gray_coding[16][i], gray_coding[16][j])
    # pprint.pprint(dist_4)
    # pprint.pprint(dist_8)

def report_ber(filename_prefix, level_list, hint=None):
    res = []
    for i in level_list:
        dist = 0
        num_bits = 0
        if i == 4:
            dist = dist_4
            num_bits = 2
        elif i == 8:
            dist = dist_8
            num_bits = 3
        elif i == 16:
            dist = dist_16
            num_bits = 4
        else:
            assert False
        fname = directory + filename_prefix + str(i)
        matrix = get_matrix_from_file(fname)
        ber_matrix = np.multiply(matrix, dist) / i
        # pprint.pprint(ber_matrix)
        ber_avg = np.sum(ber_matrix) / num_bits
        if hint is None:
            print("'" + filename_prefix + str(i) + "' :", str(ber_avg)+",")
        else:
            print("'" + hint + str(i) + "' :", str(ber_avg)+",")
        res.append(ber_avg)
    return res

def report_ber_sample(filename_prefix, level_list, sample_sizes, hint=None, sample_method="random"):
    res = []
    skip = []
    for i in level_list:
        dist = 0
        num_bits = 0
        if i == 4:
            dist = dist_4
            num_bits = 2
        elif i == 8:
            dist = dist_8
            num_bits = 3
        elif i == 16:
            dist = dist_16
            num_bits = 4
        else:
            assert False
        for sample_size in sample_sizes:
            try:
                fname = directory + filename_prefix + str(i) + "_" + str(sample_size) + "_" + sample_method
                matrix = get_matrix_from_file(fname)
                ber_matrix = np.multiply(matrix, dist) / i
                # pprint.pprint(ber_matrix)
                ber_avg = np.sum(ber_matrix) / num_bits
                if hint is None:
                    print("'" + filename_prefix + str(i) + "_" + str(sample_size) + "' :", str(ber_avg)+",")
                else:
                    print("'" + hint + str(i) + "_" + str(sample_size) + "' :", str(ber_avg)+",")
                res.append(ber_avg)
            except FileNotFoundError as e:
                # print(e)
                # print("filename", filename_prefix + str(i) + "_" + str(sample_size))
                skip.append(sample_size)
                continue
    return res, skip

def report_ber_reduction_sample(our, pba, hint, sample_sizes):
    # assert len(our) == len(pba)
    reductions = {}
    for i in range(0, len(hint)):
        sample_reduction = []
        for sample_size in sample_sizes:
            try:
                pba_ber = pba[sample_size]
                our_ber = our[sample_size]
                reduction = (pba_ber - our_ber) / pba_ber
                sample_reduction.append(reduction)
                print(f"BER reduction for level{hint[i]}_{sample_size} =", reduction)
            except KeyError as e:
                # print(e)
                continue
        reductions[hint[i]] = sample_reduction
    return reductions
        
def report_ber_reduction(our, pba, hint):
    assert len(our) == len(pba)
    for i in range(0, len(our)):
        print(f"BER reduction for level{hint[i]} =", (pba[i] - our[i]) / pba[i])

def trans_sample(level_list, sample_sizes, sample_method="random"):
    init_dist()
    print("raw_ber = {\\")
    pba_ber, pba_skip = report_ber_sample("dala", level_list, sample_sizes, sample_method)
    fpba_ber, fpba_skip = report_ber_sample("flexible", level_list, sample_sizes, sample_method)
    print("}")
    
    pba_sample_sizes = [item for item in sample_sizes if item not in pba_skip]
    fpba_sample_sizes = [item for item in sample_sizes if item not in fpba_skip]
    pba_dict = dict(zip(pba_sample_sizes, pba_ber))
    fpba_dict = dict(zip(fpba_sample_sizes, fpba_ber))

    reductions = report_ber_reduction_sample(fpba_dict, pba_dict, list(map(str, level_list)), sample_sizes)
    sample_sizes = list(set(pba_sample_sizes) & set(fpba_sample_sizes))
    for i in level_list:
        plt.plot(sample_sizes, reductions[str(i)], marker='o')
        plt.xlabel('sample sizes')
        plt.ylabel('BER reductions')
        plt.title(f'{i}-level BER reductions')
        plt.show()

def trans(level_list):
    init_dist()
    print("raw_ber = {\\")
    pba_ber = report_ber("dala", level_list)
    fpba_ber = report_ber("flexible", level_list)
    
    # norm_ber = report_ber("SBAmeanvar", [4, 8], hint="norm")
    print("}")
    report_ber_reduction(fpba_ber, pba_ber, list(map(str, level_list)))


if __name__ == "__main__":
    # ours_drift = report_results("flexible_dala", "pba_res")
    # pba_drift = report_results("dala", "our_res")
    # our_sigma = report_results("SBAvar", "pba_our_search")
    # our_norm = report_results("SBAmeanvar", "pba_our_search_mean")
    # report_drift_reduction([(ours_drift, pba_drift, "Overall"), 
    #                         (ours_drift, our_sigma, "RDR"), 
    #                         (ours_drift, our_norm, "Non-Normal")])
# we should use this file for final results reported in the paper
# instead of scheme_analyze.py (which is non-uniform weighted average)
    trans_sample([8, 16], range(15, 100), sample_method="random")
    # trans([8, 16])
