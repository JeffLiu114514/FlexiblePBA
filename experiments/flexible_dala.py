import json
import copy
import pprint
import bisect
from refine_utils import get_ber_for_allocs, init_dist

date="Feb12"
MAX_RES = 64
PRINT_ANCHOR = False
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

def find_leftmost(levels):
    sorted_levels = sorted(levels, key=lambda x: x[1])
    return sorted_levels[0], sorted_levels[0][1]
    
def update(R, anchor, xl, xh, BER):
    left_perc = bisect.bisect_left(R, anchor) / len(R)
    if left_perc > BER:
        return xl, 2*MAX_RES
    num_discard = int((BER - left_perc) * len(R))
    if num_discard == 0:
        return anchor, R[-1] + 1
    return anchor, R[-num_discard] + 1
    

def minimal_BER(specified_levels, eps, distributions, low_BER = 0, high_BER = 1, flexible_refine_flag=False, double=False):
    # rationale for double: for 4 levels with insufficient data to characterize the error
    #   we need to allocate 8 levels then half the levels
    if double:
        specified_levels = specified_levels * 2
    
    # ratios_lists = []
    while high_BER - low_BER > eps:
        cur_BER = (low_BER + high_BER) / 2
        
        # flexible greedy algorithm
        result_level = []
        cur_levels = candidate_gen(cur_BER, distributions)
        confirmed_level, anchor = find_leftmost(cur_levels)
        if PRINT_ANCHOR:
            print("confirmed level and anchor:", confirmed_level[0], confirmed_level[1], anchor)
        result_level.append([confirmed_level[0], confirmed_level[1], confirmed_level[3], confirmed_level[4]])
        while anchor <= MAX_RES+1:
            temp_levels = []
            for cur_level in cur_levels:
                Rlow, Rhigh, RelaxDistr, tmin, tmax = cur_level
                if Rlow < anchor:
                    Rlow, Rhigh = update(RelaxDistr, anchor, Rlow, Rhigh, cur_BER)
                if Rlow >= anchor:
                    temp_levels.append([Rlow, Rhigh, RelaxDistr, tmin, tmax])
            if len(temp_levels) == 0:
                break
            confirmed_level, anchor = find_leftmost(temp_levels)
            # anchor += 1
            
            if PRINT_ANCHOR:
                print("confirmed level and anchor:", confirmed_level[0], confirmed_level[1],confirmed_level[3],confirmed_level[4], anchor)

            result_level.append([confirmed_level[0], confirmed_level[1], confirmed_level[3], confirmed_level[4]])
            cur_levels = temp_levels   
            
        
        cur_levels = result_level
        if DEBUG: print(len(cur_levels), cur_BER)
        if len(cur_levels) < specified_levels: # the precision requirement is too strict to be met
            low_BER = cur_BER # make next BER bigger
        elif len(cur_levels) > specified_levels:
            high_BER = cur_BER
        else:
            high_BER = cur_BER
            best_level, best_BER = cur_levels, cur_BER
            
            
            # if flexible_refine_flag:
            #     # print("current gamma: ", cur_BER)
            #     # print("current level: ", cur_levels)
            #     refined, min_ber = flexible_refine(cur_levels, specified_levels, distributions)
            #     ratios_lists.append({"gamma": cur_BER, "level": refined, "min_ber": min_ber})
            # else:
            #     refined = refine(cur_levels)
        
    if double:
        best_level = half(best_level)
    
    if flexible_refine_flag:
        # print("current gamma: ", cur_BER)
        refined = flexible_refine(best_level, specified_levels, distributions)
    else:
        refined = refine(best_level)

    # for item in enumerate(ratios_lists):
    #     print(item)
    #     print()
    
    if DEBUG: print(refined, best_BER)
    assert len(refined) == specified_levels / 2 if double else specified_levels
    return refined, best_BER

def candidate_gen(BER, distributions):
    if DEBUG:
        print(BER, "Started")
    levels = []
    flag = False
    for tmin in range(0, 60):
        tmax = tmin + 4
        RelaxDistr = distributions[(tmin, tmax)]
        # if DEBUG:
        #     print(len(RelaxDistr),int(BER * len(RelaxDistr) / 2))
        if int(BER * len(RelaxDistr)) == 0: flag = True
        Rlow, Rhigh = getReadRange(RelaxDistr, BER)
        # assert Rlow <= tmin and tmax <= Rhigh, (Rlow, Rhigh, tmin, tmax)
        levels.append([Rlow, Rhigh, RelaxDistr, tmin, tmax])
    if DEBUG: print("BER * length of distribution is 0: ", flag)
    return levels


def getReadRange(vals, BER):
    # The read range [Rmin, Rmax) -- any point within [Rmin, Rmax) are within this level
    num_discard = int(BER * len(vals))
    if num_discard == 0:
        return vals[0], vals[-1] + 1
    return vals[0], vals[-num_discard] + 1


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
    print("vanilla: ", vanilla)
    print("naive refine: ", level_alloc)
    
    dist_4, dist_8, dist_16 = init_dist()
    
    for i in range(1, len(vanilla)):
        assert level_alloc[i - 1][1] <= level_alloc[i][0]
        min_ber = get_ber_for_allocs(level_alloc, distributions, specified_levels, dist_4, dist_8, dist_16)
        best_j = level_alloc[i - 1][1]
        for j in range(vanilla[i - 1][0], vanilla[i][1]):
            level_alloc[i - 1][1] = j
            level_alloc[i][0] = j
            ber = get_ber_for_allocs(level_alloc, distributions, specified_levels, dist_4, dist_8, dist_16)
            # print(j, ber)
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

    # print("flexible refine", level_alloc)
    # print("BER: ", min_ber)
    return level_alloc, min_ber

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
    write_to_json(bpc, f"tests/{bits_per_cell}bpc_flexible_dala_{date}.json")


if __name__ == "__main__":
    distributions = init_model()
    # refined, best_BER = minimal_BER(4, 1e-3, distributions, 0, 1, True, True)
    # dump_to_json(refined)
    # refined, best_BER = minimal_BER(8, 1e-3, distributions, flexible_refine_flag=False)
    # print(refined, best_BER)
    refined, best_BER = minimal_BER(8, 1e-10, distributions, flexible_refine_flag=True)
    print(refined, best_BER)
    
    # refined, best_BER = minimal_BER(16, 1e-3, distributions, flexible_refine_flag=False)
    # print(refined, best_BER)
    # refined, best_BER = minimal_BER(16, 1e-3, distributions, flexible_refine_flag=True)
    # print(refined, best_BER)
    
    # dump_to_json(refined)
    # dump_to_json(minimal_BER(16, 1e-10))
