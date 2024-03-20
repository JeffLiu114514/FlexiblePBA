import dala
import flexible_dala
import numpy

outfile = "./ember_capacity/"

def decide_end_level(point, level_alloc):
    for i in range(len(level_alloc)):
        rmin, rmax, wmin, wmax = level_alloc[i]    
        if point >= rmin and point < rmax:
            return i
    if point == 64: # last level
        return len(level_alloc)-1
    assert False, point

def simulate_error(level_alloc, distributions):
    num_levels = len(level_alloc)
    P = numpy.zeros((num_levels, num_levels))
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

def get_dala(distributions):
    res = {}
    for i in [4,8,16]:
        if i <= 5:
            res[i] = dala.minimal_BER(i, 1e-3, distributions, 0, 1, True)
        else:
            res[i] = dala.minimal_BER(i, 1e-3, distributions)
    # print(res)
    return res

def get_flexible_dala(distributions):
    res = {}
    for i in [4,8,16]:
        if i <= 5:
            res[i] = flexible_dala.minimal_BER(i, 1e-3, distributions, 0, 1, True)
        else:
            res[i] = flexible_dala.minimal_BER(i, 1e-3, distributions)
    # print(res)
    return res

def simulate_all_levels(dala_allocs, distribution, outfile):
    for i in [4,8,16]:
        P = simulate_error(dala_allocs[i][0], distribution)
        dump_matrix(P, outfile)

def dump_matrix(matrix, hint):
    num_level = len(matrix)
    with open(hint + str(num_level), "w") as fout:
        to_write = []
        for i in range(len(matrix)):
            to_write.append(",".join(map(str, matrix[i])) + "\n")
        fout.writelines(to_write)
        
def dump_matrix_sample(matrix, hint, sample_size, suffix=""):
    num_level = len(matrix)
    with open(hint + str(num_level) + "_" + str(sample_size) + "_" + suffix, "w") as fout:
        to_write = []
        for i in range(len(matrix)):
            to_write.append(",".join(map(str, matrix[i])) + "\n")
        fout.writelines(to_write)

if __name__ == "__main__":
    distributions = dala.init_model()
    dala_allocs = get_dala(distributions)
    # print(dala_allocs[16])
    fdala_allocs = get_flexible_dala(distributions)
    # print(fdala_allocs[16])
    simulate_all_levels(dala_allocs, distributions, outfile+"dala")
    simulate_all_levels(fdala_allocs, distributions, outfile+"flexible")
