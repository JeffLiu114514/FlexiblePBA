from dala import minimal_BER as dala_minimal_BER
from dala import getReadRange, refine, flexible_refine
from flexible_dala import minimal_BER as flexible_dala_minimal_BER
from dala_genmatrix import simulate_error
from ecc import bestcode, allcode
import numpy as np
import networkx as nx
import numpy as np
import copy
import json
from tqdm import tqdm
import argparse

DEBUG=False
ecc_params = {"codes": allcode(),
              "spec_ber": 1e-14,
              "maxk": 2**12, 
              "maxn": 2**12}

# Basic functions
def init_model(model_name="retention1s.csv"):
    distributions = {}
    with open("../model/"+model_name, "r") as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.split(',')
            tmin, tmax, distr = int(tokens[0]), int(tokens[1]), list(map(int, tokens[2:]))
            distributions[(tmin, tmax)] = distr
    # print(distributions)
    return distributions

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

def get_ber_for_allocs(dala_alloc, distributions, n):
    dist_4, dist_8, dist_16 = init_dist()
    P = simulate_error(dala_alloc, distributions)
    ber = report_ber(P, n, dist_4, dist_8, dist_16)
    return ber

def get_ecc_for_ber(ber):
    ecc = bestcode(ecc_params["codes"], ecc_params["spec_ber"], ber, ecc_params["maxn"], ecc_params["maxk"])
    return ecc[1]



# Graph based functions
def dala_level_inference(BER, distributions):
    if DEBUG:
        print(BER, "Started")
    levels = []
    flag = False
    for tmin in range(0, 60):
        tmax = tmin + 4
        RelaxDistr = distributions[(tmin, tmax)]
        if DEBUG:
            print(len(RelaxDistr),int(BER * len(RelaxDistr) / 2))
        if DEBUG: 
            if int(BER * len(RelaxDistr) / 2) == 0: flag = True
        Rlow, Rhigh = getReadRange(RelaxDistr, BER)
        # assert Rlow <= tmin and tmax <= Rhigh, (Rlow, Rhigh, tmin, tmax)
        levels.append([Rlow, Rhigh, tmin, tmax])
    if DEBUG: print("BER * length of distribution/2 is 0: ", flag)
    # print(levels)
    return np.unique(np.array(levels), axis=0).tolist()

def fdala_level_inference(BER, distributions):
    levels = []
    for tmin in range(0, 60):
        tmax = tmin + 4
        RelaxDistr = distributions[(tmin, tmax)]
        num_discard = int(BER * len(RelaxDistr))
        if num_discard == 0:
            levels.append([RelaxDistr[0], RelaxDistr[-1] + 1, tmin, tmax])
        else:
            # print("num_discard: ", num_discard)
            for i in range(0, num_discard):
                Rlow, Rhigh = RelaxDistr[i], RelaxDistr[-(num_discard-i)]+1
                levels.append([Rlow, Rhigh, tmin, tmax])
    # print(levels)
    return np.unique(np.array(levels), axis=0).tolist()

def leftmost_readrange(vals, BER):
    # The read range [Rmin, Rmax) -- any point within [Rmin, Rmax) are within this level
    num_discard = int(BER * len(vals))
    if num_discard == 0:
        return vals[0], vals[-1] + 1
    return vals[0], vals[-num_discard] + 1

def is_overlap(a, b):
    overlap = max(0, min(a[1], b[1]) - max(a[0], b[0]))
    if overlap == 0:
        return False
    return True

def construct_graph(levels):
    graph = nx.Graph()
    
    # Add levels as nodes with their index as id
    for i, level in enumerate(levels):
        graph.add_node(i, level=level)
    
    # Iterate over all combinations of any pair of levels
    id_to_levels = nx.get_node_attributes(graph, "level")
    for i in range(len(levels)):
        for j in range(i+1, len(levels)):
            level1 = id_to_levels[i]
            level2 = id_to_levels[j]
            
            # Check if levels overlap
            if not is_overlap(level1, level2):
                graph.add_edge(i, j)
                
    return graph

def dala_graph(specified_levels, cur_gamma, distributions, flexible_refine_flag=False):
    levels = dala_level_inference(cur_gamma, distributions)
    graph = construct_graph(levels)
    num_nodes = graph.number_of_nodes()
    print("Number of nodes in the graph:", num_nodes)
    # cliques = list(nx.enumerate_all_cliques(graph))
    cliques = list(nx.find_cliques(graph))
    # print("num cliques: ", len(cliques))
    cliques_ = [clique for clique in cliques if len(clique) == specified_levels]
    # print("num cliques with", specified_levels, "levels:", len(cliques_))
    
    best_ber = 1
    best_clique = None
    id_to_levels = nx.get_node_attributes(graph, "level")
    for clique in tqdm(cliques_):
        # print("clique: ", clique)
        level_alloc = []
        for i in clique:
            level = id_to_levels[i]
            # print(level)
            level_alloc.append(copy.deepcopy(level))
        level_alloc.sort(key=lambda x: x[0])
        if flexible_refine_flag:
            refined = flexible_refine(level_alloc, specified_levels, distributions)
        else:
            refined = refine(level_alloc)
        ber = get_ber_for_allocs(refined, distributions, specified_levels)
        if ber < best_ber:
            best_ber = ber
            best_clique = refined
    
    # print(best_clique)
    # print(best_ber)
    return best_clique, best_ber

def fdala_graph(specified_levels, cur_gamma, distributions, flexible_refine_flag=False):

    levels = fdala_level_inference(cur_gamma, distributions)
    graph = construct_graph(levels)
    num_nodes = graph.number_of_nodes()
    # print("Number of nodes in the graph:", num_nodes)
    # cliques = list(nx.enumerate_all_cliques(graph))
    cliques = list(nx.find_cliques(graph))
    # print("num cliques: ", len(cliques))
    cliques_ = [clique for clique in cliques if len(clique) == specified_levels]
    # print("num cliques with", specified_levels, "levels:", len(cliques_))
    
    best_ber = 1
    best_clique = None
    id_to_levels = nx.get_node_attributes(graph, "level")
    
    for clique in tqdm(cliques_):
        # print("clique: ", clique)
        level_alloc = []
        for i in clique:
            level = id_to_levels[i]
            # print(level)
            level_alloc.append(copy.deepcopy(level))
        level_alloc.sort(key=lambda x: x[0])
        if flexible_refine_flag:
            refined = flexible_refine(level_alloc, specified_levels, distributions)
        else:
            refined = refine(level_alloc)
        ber = get_ber_for_allocs(refined, distributions, specified_levels)
        if ber < best_ber:
            best_ber = ber
            best_clique = refined
    
    # print(best_clique)
    # print(best_ber)
    return best_clique, best_ber


# Test functions
def run_basic_test(distributions, num_levels, eps):
    print("Basic tests for", num_levels, "levels")
    n = num_levels
    results = {}

    dala_refined, gamma = dala_minimal_BER(n, eps, distributions)
    ber = get_ber_for_allocs(dala_refined, distributions, n)
    ecc = get_ecc_for_ber(ber)
    results["dala"] = (ber, ecc, gamma)
    
    dala_refined, gamma = dala_minimal_BER(n, eps, distributions, flexible_refine_flag=True)
    ber = get_ber_for_allocs(dala_refined, distributions, n)
    results["dala+flexible_refine"] = (ber, ecc, gamma)
    
    flexible_refined, gamma, _ = flexible_dala_minimal_BER(n, eps, distributions)
    ber = get_ber_for_allocs(flexible_refined, distributions, n)
    results["flexible_dala"] = (ber, ecc, gamma)
    
    flexible_refined, gamma, _ = flexible_dala_minimal_BER(n, eps, distributions, flexible_refine_flag=True)
    ber = get_ber_for_allocs(flexible_refined[0], distributions, n)
    results["flexible_dala+flexible_refine"] = (ber, ecc, gamma)

    return results

def run_graph_test(distributions, num_levels, dala_gamma, fdala_gamma):
    print("Graph tests for", num_levels, "levels")
    n = num_levels
    results = {}
    
    dala_clique, dala_ber = dala_graph(n, dala_gamma, distributions)
    ecc = get_ecc_for_ber(dala_ber)
    results["dala_graph"] = (dala_ber, ecc)
    
    if n <= 8:
        fdala_clique, fdala_ber = fdala_graph(n, fdala_gamma, distributions)
        ecc = get_ecc_for_ber(fdala_ber)
        results["fdala_graph"] = (fdala_ber, ecc)
    
    return results

def run_graph_relaxation_test(distributions, num_levels, dala_gamma, fdala_flag=False, fdala_gamma=None):
    print("Graph relaxation tests for", num_levels, "levels")
    results = {}
    min_ber = {"ber": 1, "ecc": 100, "clique": None, "relaxation": None}
    
    for relaxation in range(11, 31, 1):
        gamma_relaxation = relaxation / 10
        best_clique, best_ber = dala_graph(num_levels, dala_gamma*gamma_relaxation, distributions)
        ecc = get_ecc_for_ber(best_ber)
        results[f"dala_graph_relaxation_{gamma_relaxation}"] = (best_ber, ecc)
        if best_ber < min_ber["ber"]:
            min_ber["ber"] = best_ber
            min_ber["ecc"] = ecc
            min_ber["clique"] = best_clique
            min_ber["relaxation"] = gamma_relaxation
    results["dala_graph_relaxation_best"] = (min_ber["ber"], min_ber["ecc"], min_ber["relaxation"])
      
    if fdala_flag:
        for relaxation in range(11, 31, 1):
            gamma_relaxation = relaxation / 10
            best_clique, best_ber = fdala_graph(num_levels, fdala_gamma*gamma_relaxation, distributions)
            ecc = get_ecc_for_ber(best_ber)
            results[f"fdala_graph_relaxation_{gamma_relaxation}"] = (best_ber, ecc)
            if best_ber < min_ber["ber"]:
                min_ber["ber"] = best_ber
                min_ber["ecc"] = ecc
                min_ber["clique"] = best_clique
                min_ber["relaxation"] = gamma_relaxation
        results["fdala_graph_relaxation_best"] = (min_ber["ber"], min_ber["ecc"], min_ber["relaxation"])
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--model_filename", type=str, help="filename for the test model file", default="retention1s.csv")
    # parser.add_argument("-n", "--num_levels", type=int, help="number of allocation levels", default=8)
    parser.add_argument("-e", "--eps", type=float, help="minimum granularity for gamma", default=1e-3)
    args = parser.parse_args()
    model_filename = args.model_filename
    # num_levels = args.num_levels
    eps = args.eps
    spec_ber = 1e-14
    
    distributions = init_model(model_filename)
    
    
    results8 = {}
    
    basic = run_basic_test(distributions, 8, eps)
    results8.update(basic)
    
    graph = run_graph_test(distributions, 8, basic["dala"][2], basic["flexible_dala"][2])
    results8.update(graph)
    
    # relaxation = run_graph_relaxation_test(distributions, 8, basic["dala"][2], True, basic["flexible_dala"][2])
    # results8.update(relaxation)
    
    with open("./all_tests/"+model_filename+"_8levels.json", "w") as f:
        json.dump(results8, f)
        
    results16 = {}
    
    basic = run_basic_test(distributions, 16, eps)
    results16.update(basic)
    
    graph = run_graph_test(distributions, 16, basic["dala"][2], basic["flexible_dala"][2])
    results16.update(graph)
    
    relaxation = run_graph_relaxation_test(distributions, 16, basic["dala"][2], False)
    results16.update(relaxation)
    
    with open("./all_tests/"+model_filename+"_16levels.json", "w") as f:
        json.dump(results16, f)
    
    