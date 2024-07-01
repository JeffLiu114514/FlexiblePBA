from dala import init_model, getReadRange, refine, flexible_refine, minimal_BER
from flexible_dala import minimal_BER as flexible_dala_minimal_BER
import networkx as nx
import numpy as np
import copy
from allinone import get_ber_for_allocs
from tqdm import tqdm

DEBUG = False

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


def dala_graph(specified_levels, cur_gamma, distributions):
    # cur_gamma = (low_gamma + high_gamma) / 2
    # while high_gamma - low_gamma > eps:
    #     print("gamma: ", cur_gamma)
    #     levels = dala_level_inference(cur_gamma, distributions)
    #     graph = construct_graph(levels)
    #     cliques = list(nx.enumerate_all_cliques(graph))
    #     if any(len(clique) >= specified_levels for clique in cliques):
    #         cur_gamma = (low_gamma + cur_gamma) / 2
    #         continue
    #     elif all(len(clique) < specified_levels for clique in cliques):
    #         cur_gamma = (cur_gamma + high_gamma) / 2
    #         continue
        
    #     print("gamma: ", cur_gamma)
    #     break
    levels = dala_level_inference(cur_gamma, distributions)
    graph = construct_graph(levels)
    num_nodes = graph.number_of_nodes()
    print("Number of nodes in the graph:", num_nodes)
    # cliques = list(nx.enumerate_all_cliques(graph))
    cliques = list(nx.find_cliques(graph))
    print("num cliques: ", len(cliques))
    cliques_ = [clique for clique in cliques if len(clique) == specified_levels]
    print("num cliques with", specified_levels, "levels:", len(cliques_))
    
    best_ber = 1
    best_clique = None
    id_to_levels = nx.get_node_attributes(graph, "level")
    for clique in tqdm(cliques_):
        # print("clique: ", clique)
        level_alloc = []
        for i in clique:
            level = id_to_levels[i]
            # print(level)
            # tmin = i
            # tmax = tmin + 4
            # RelaxDistr = distributions[(tmin, tmax)]
            # if DEBUG:
            #     print(len(RelaxDistr),int(cur_gamma * len(RelaxDistr) / 2))
            # # if DEBUG: 
            # #     if int(cur_gamma * len(RelaxDistr) / 2) == 0: flag = True
            # Rlow, Rhigh = getReadRange(RelaxDistr, cur_gamma)
            # level_alloc.append([Rlow, Rhigh, tmin, tmax])
            level_alloc.append(copy.deepcopy(level))
        level_alloc.sort(key=lambda x: x[0])
        refined = refine(level_alloc) # refine(level_alloc)  flexible_refine(level_alloc, specified_levels, distributions)
        ber = get_ber_for_allocs(refined, distributions, specified_levels)
        if ber < best_ber:
            best_ber = ber
            best_clique = refined
    
    print(best_clique)
    print(best_ber)
    return best_clique, best_ber

def fdala_graph(specified_levels, cur_gamma, distributions):

    levels = fdala_level_inference(cur_gamma, distributions)
    graph = construct_graph(levels)
    num_nodes = graph.number_of_nodes()
    print("Number of nodes in the graph:", num_nodes)
    # cliques = list(nx.enumerate_all_cliques(graph))
    cliques = list(nx.find_cliques(graph))
    print("num cliques: ", len(cliques))
    cliques_ = [clique for clique in cliques if len(clique) == specified_levels]
    print("num cliques with", specified_levels, "levels:", len(cliques_))
    
    best_ber = 1
    best_clique = None
    id_to_levels = nx.get_node_attributes(graph, "level")
    
    for clique in tqdm(cliques_):
        # print("clique: ", clique)
        level_alloc = []
        for i in clique:
            level = id_to_levels[i]
            # print(level)
            # tmin = i
            # tmax = tmin + 4
            # RelaxDistr = distributions[(tmin, tmax)]
            # if DEBUG:
            #     print(len(RelaxDistr),int(cur_gamma * len(RelaxDistr) / 2))
            # # if DEBUG: 
            # #     if int(cur_gamma * len(RelaxDistr) / 2) == 0: flag = True
            # Rlow, Rhigh = getReadRange(RelaxDistr, cur_gamma)
            # level_alloc.append([Rlow, Rhigh, tmin, tmax])
            level_alloc.append(copy.deepcopy(level))
        level_alloc.sort(key=lambda x: x[0])
        refined = flexible_refine(level_alloc, specified_levels, distributions)#refine(level_alloc)
        ber = get_ber_for_allocs(refined, distributions, specified_levels)
        if ber < best_ber:
            best_ber = ber
            best_clique = refined
    
    print(best_clique)
    print(best_ber)
    return best_clique, best_ber



if __name__ == "__main__":
    distributions = init_model()
    # gamma_relaxation = 1.5
    
    num_levels = 16
    
    refined, dala_best_gamma = minimal_BER(num_levels, 1e-3, distributions)
    ber = get_ber_for_allocs(refined, distributions, num_levels)
    print("original dala")
    print("dala alloc: ", refined)
    print("dala ber: ", ber)
    
    refined, fdala_best_gamma, _ = flexible_dala_minimal_BER(num_levels, 1e-3, distributions)
    ber = get_ber_for_allocs(refined, distributions, num_levels)
    print("flexible dala")
    print("fdala alloc: ", refined)
    print("dala ber: ", ber)
    
    print("graph based dala")
    min_ber = {"ber": 1, "clique": None, "relaxation": None}
    best_clique, best_ber = dala_graph(num_levels, dala_best_gamma, distributions)
    # for relaxation in range(11, 31, 1):
    #     gamma_relaxation = relaxation / 10
    #     print("relaxation: ", gamma_relaxation)
    #     best_clique, best_ber = dala_graph(8, dala_best_gamma*gamma_relaxation, distributions)
    #     if best_ber < min_ber["ber"]:
    #         min_ber["ber"] = best_ber
    #         min_ber["clique"] = best_clique
    #         min_ber["relaxation"] = gamma_relaxation
            
    # print("graph based flexible dala")
    # min_ber = {"ber": 1, "clique": None, "relaxation": None}
    # best_clique, best_ber = fdala_graph(num_levels, fdala_best_gamma, distributions)
    # for relaxation in range(11, 12, 1):
    #     gamma_relaxation = relaxation / 10
    #     print("relaxation: ", gamma_relaxation)
    #     best_clique, best_ber = fdala_graph(8, best_gamma*gamma_relaxation, distributions)
    #     if best_ber < min_ber["ber"]:
    #         min_ber["ber"] = best_ber
    #         min_ber["clique"] = best_clique
    #         min_ber["relaxation"] = gamma_relaxation
    
    # print("best clique: ", min_ber["clique"])
    # print("best ber: ", min_ber["ber"])
    # print("best relaxation: ", min_ber["relaxation"])
    
    # best_clique, best_ber = minimal_BER_graph(8, best_gamma*gamma_relaxation, distributions)
