import json
import copy
import pprint
from refine_utils import get_ber_for_allocs, init_dist

date="Feb6"

DEBUG = False


def init_model():
    distributions = {}
    with open("../model/retention1s.csv", "r") as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.split(',')
            tmin, tmax, distr = int(tokens[0]), int(tokens[1]), list(map(int, tokens[2:]))
            distributions[(tmin, tmax)] = distr
    # print(distributions)
    return distributions

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
    return longest_non_overlap(levels)

def longest_non_overlap(levels):
    # this is a greedy algorithm
    # levels is a list of levels, each level is a list of [Rlow, Rhigh, tmin, tmax]
    res = []
    # first sort by Rhigh
    sorted_levels = sorted(levels, key=lambda x: x[1])
    res.append(sorted_levels[0])
    cur = sorted_levels[0]
    for i in range(1, len(sorted_levels)):
        nxt = sorted_levels[i]
        # the next level's Rlow does not overlap with the current level's Rhigh
        # the next level's tmin (write ranges) does not overlap with current level's tmax
        if nxt[0] >= cur[1]:# and nxt[2] >= cur[3]
            res.append(nxt)
            cur = nxt
    return res


def getReadRange(vals, BER):
    # The read range [Rmin, Rmax) -- any point within [Rmin, Rmax) are within this level
    num_discard = int(BER * len(vals) / 2)
    if num_discard == 0:
        return vals[num_discard], vals[-1] + 1
    return vals[num_discard], vals[-num_discard] + 1
    # total_num_discard = int(BER * len(vals))
    # if total_num_discard % 2 == 0:
    #     num_discard_front, num_discard_back = int(total_num_discard / 2), int(total_num_discard / 2)
    # else:
    #     num_discard_front, num_discard_back = int(total_num_discard / 2) + 1, int(total_num_discard / 2)
    # return vals[num_discard_front], vals[-num_discard_back-1] + 1

def refine(level_alloc):
    '''
    close the gap between adjacent read ranges
    Example: list of [Rlow, Rhigh, tmin, tmax]
        [2, 14, 0, 4], [16, 28, 16, 20], [32, 44, 33, 37], [48, 56, 46, 50]
    --> [2, 15, 0, 4], [15, 30, 16, 20], [30, 46, 33, 37], [46, 56, 46, 50]
    --> [0, ...                                               , 63, ...  ]
    '''
    # print("before refine", level_alloc)
    for i in range(1, len(level_alloc)):
        assert level_alloc[i - 1][1] <= level_alloc[i][0] 
        merge = int((level_alloc[i - 1][1] + level_alloc[i][0]) / 2)
        level_alloc[i - 1][1] = merge
        level_alloc[i][0] = merge
    level_alloc[0][0] = 0
    level_alloc[len(level_alloc)-1][1] = 64
    return level_alloc

def flexible_refine(level_alloc, specified_levels, distributions):
    '''
    optimally close the gap between adjacent read ranges with respect to BER
    '''
    
    vanilla = copy.deepcopy(level_alloc)
    level_alloc = refine(level_alloc)
    if DEBUG:
        print("vanilla: ", vanilla)
        print("naive refine: ", level_alloc)
    
    dist_4, dist_8, dist_16 = init_dist()
    
    for i in range(1, len(vanilla)):
        assert level_alloc[i - 1][1] <= level_alloc[i][0]
        min_ber = get_ber_for_allocs(level_alloc, distributions, specified_levels, dist_4, dist_8, dist_16)
        best_j = level_alloc[i - 1][1]
        for j in range(vanilla[i - 1][0], vanilla[i][1]+1):
            level_alloc[i - 1][1] = j
            level_alloc[i][0] = j
            ber = get_ber_for_allocs(level_alloc, distributions, specified_levels, dist_4, dist_8, dist_16)
            if DEBUG:print(j, ber)
            # print(level_alloc)
            if ber < min_ber:
                min_ber = ber
                best_j = j
            else:
                level_alloc[i - 1][1] = best_j
                level_alloc[i][0] = best_j
                continue
    level_alloc[0][0] = 0
    level_alloc[len(level_alloc)-1][1] = 64

    if DEBUG:
        print("flexible refine", level_alloc)
        print("BER: ", min_ber)
    return level_alloc


def half(level_alloc):
    assert len(level_alloc) % 2 == 0
    res = []
    '''
    [Rlow_i, Rhigh_i, tmin_i, tmax_i], [Rlow_{i+1}, Rhigh_{i+1}, tmin_{i+1}, tmax_{i+1}]
    -->
    [Rlow_i, Rhigh_{i+1}, (tmin_i+tmin_{i+1} ) / 2, (tmax_i+tmax_{i+1} ) / 2] 
    '''
    for i in range(0, int(len(level_alloc) / 2)):
        res.append([level_alloc[i*2][0], level_alloc[i*2+1][1],
                    int((level_alloc[i*2][2] + level_alloc[i*2+1][2]) / 2),
                    int((level_alloc[i*2][3] + level_alloc[i*2+1][3]) / 2)])
    assert len(res) == len(level_alloc) / 2
    return res



def minimal_BER(specified_levels, eps, distributions, low_BER = 0, high_BER = 1, flexible_refine_flag=False, double=False):
    # rationale for double: for 4 levels with insufficient data to characterize the error
    #   we need to allocate 8 levels then half the levels
    if double:
        specified_levels = specified_levels * 2
    while high_BER - low_BER > eps:
        cur_BER = (low_BER + high_BER) / 2
        cur_levels = level_inference(cur_BER, distributions)
        if DEBUG: print(len(cur_levels), cur_BER)
        if len(cur_levels) < specified_levels: # the precision requirement is too strict to be met
            low_BER = cur_BER # make next BER bigger
        elif len(cur_levels) > specified_levels:
            high_BER = cur_BER
        else:
            high_BER = cur_BER
            best_level, best_BER = cur_levels, cur_BER
    if double:
        best_level = half(best_level)
    if flexible_refine_flag:
        refined = flexible_refine(best_level, specified_levels, distributions)
    else:
        refined = refine(best_level)
    if DEBUG: print(refined, best_BER)
    assert len(refined) == specified_levels / 2 if double else specified_levels
    return refined, best_BER

def read_from_json(filename):
    return json.load(open(filename))

def write_to_json(data, filename):
    json.dump(data, open(filename, "w"), indent=1)

def dump_to_json(level_alloc):
    if len(level_alloc) == 16:
        bits_per_cell = 4
    elif len(level_alloc) == 8:
        bits_per_cell = 3
    elif len(level_alloc) == 4:
        bits_per_cell = 2
    bpc = read_from_json(f"settings/{bits_per_cell}bpc.json")
    for i in range(0, len(level_alloc)):
        # [Rlow, Rhigh, tmin, tmax]
        bpc['level_settings'][i]["adc_upper_read_ref_lvl"] = level_alloc[i][1]
        bpc['level_settings'][i]["adc_lower_write_ref_lvl"] = level_alloc[i][2]
        bpc['level_settings'][i]["adc_upper_write_ref_lvl"] = level_alloc[i][3]
    write_to_json(bpc, f"tests/{bits_per_cell}bpc_dala_{date}.json")


if __name__ == "__main__":
    distributions = init_model()
    # refined, best_BER = minimal_BER(4, 1e-3, distributions, 0, 1, True)
    # dump_to_json(refined)
    # refined, best_BER = minimal_BER(8, 1e-3, distributions)
    # dump_to_json(refined)
    # dump_to_json(minimal_BER(16, 1e-10))
    
    refined, best_BER = minimal_BER(8, 1e-3, distributions, flexible_refine_flag=False)
    print(refined, best_BER)
    # refined, best_BER = minimal_BER(8, 1e-3, distributions, flexible_refine_flag=True)
    # print(refined, best_BER)
    
    # refined, best_BER = minimal_BER(16, 1e-3, distributions, flexible_refine_flag=False)
    # print(refined, best_BER)
    # refined, best_BER = minimal_BER(16, 1e-3, distributions, flexible_refine_flag=True)
    # print(refined, best_BER)
