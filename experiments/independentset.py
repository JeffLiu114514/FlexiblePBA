from dala import init_model, getReadRange, refine, flexible_refine, minimal_BER
import networkx as nx
from allinone import get_ber_for_allocs

DEBUG = False

def level_inference(BER, distributions):
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
    return levels

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
    for i in range(len(levels)):
        for j in range(i+1, len(levels)):
            level1 = levels[i]
            level2 = levels[j]
            
            # Check if levels overlap
            if not is_overlap(level1, level2):
                graph.add_edge(i, j)
    
    return graph


def minimal_BER_graph(specified_levels, cur_gamma, distributions):
    # cur_gamma = (low_gamma + high_gamma) / 2
    # while high_gamma - low_gamma > eps:
    #     print("gamma: ", cur_gamma)
    #     levels = level_inference(cur_gamma, distributions)
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
    levels = level_inference(cur_gamma, distributions)
    graph = construct_graph(levels)
    cliques = list(nx.enumerate_all_cliques(graph))
    cliques_ = [clique for clique in cliques if len(clique) == specified_levels]
    # print(cliques_)
    
    best_ber = 1
    best_clique = None
    for clique in cliques_:
        level_alloc = []
        for i in clique:
            tmin = i
            tmax = tmin + 4
            RelaxDistr = distributions[(tmin, tmax)]
            if DEBUG:
                print(len(RelaxDistr),int(cur_gamma * len(RelaxDistr) / 2))
            if DEBUG: 
                if int(cur_gamma * len(RelaxDistr) / 2) == 0: flag = True
            Rlow, Rhigh = getReadRange(RelaxDistr, cur_gamma)
            level_alloc.append([Rlow, Rhigh, tmin, tmax])
        refined = refine(level_alloc)
        ber = get_ber_for_allocs(refined, distributions, 8)
        if ber < best_ber:
            best_ber = ber
            best_clique = refined
    
    print(best_clique)
    print(best_ber)
    return best_clique, best_ber


if __name__ == "__main__":
    distributions = init_model()
    
    refined, best_gamma = minimal_BER(8, 1e-3, distributions)
    ber = get_ber_for_allocs(refined, distributions, 8)
    print("dala ber: ", ber)
    print("dala alloc: ", refined)
    
    best_clique, best_ber = minimal_BER_graph(8, best_gamma*1.5, distributions)
    
    
    
    
    