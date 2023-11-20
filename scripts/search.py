from copy import deepcopy
from itertools import combinations
from itertools import permutations


def tree_original(g, h, i_g, i_h):
    if i_h == len(h):
        print()
        return
    if i_g == len(g):
        print()
    for i in range(i_h, len(h)):
        print(f"{g[i_g]}-{h[i]}")
        # Swap
        g_new = deepcopy(g)
        h_new = deepcopy(h)
        temp = h_new[i_h]
        h_new[i_h] = h_new[i]
        h_new[i] = temp
        tree_original(g_new, h_new, i_g + 1, i_h + 1)


def tree_pruned(g, h, i_g, i_h, D, B):
    if i_h == len(h):
        # if i_g < len(g):
        #     print(f"{g[i_g]}-{h[i_h - 1]}")
        B += 1
        print()
        return
    if i_g == len(g):
        return
    for i in range(i_h, len(h)):
        print(f"{g[i_g]}-{h[i]}")
        # Swap
        g_new = deepcopy(g)
        h_new = deepcopy(h)
        temp = h_new[i_h]
        h_new[i_h] = h_new[i]
        h_new[i] = temp
        tree_pruned(g_new, h_new, i_g + 1, i_h+1, D+1, B)


def print_list(lst):
    for i in lst:
        print(i, end="")
    print()


mcs_count = 0
new_count = 0

def mcsp_tree(g, h, i_g, i_h):
    global mcs_count
    if i_h == len(h):
        mcs_count += 1
        # print_list(g)
        # print_list(h)
        # print()
        return
    if i_g == len(g):
        print()
    for i in range(i_h, len(h)):
        # print(f"{g[i_g]}-{h[i]}")
        # Swap
        h_new = [h[i]] + h[0:i] + h[i+1:len(h)]
        mcsp_tree(g, h_new, i_g + 1, i_h + 1)


def new_tree(g, h, i_g, i_h):
    global new_count
    if i_h == len(h):
        new_count += 1
        # print_list(g)
        # print_list(h)
        # print()
        return
    if i_g == len(g):
        print()
    for i in range(i_h, len(h)):
        # print(f"{g[i_g]}-{h[i]}")
        # Swap
        h_new = [h[i]] + h[0:i] + h[i+1:len(h)]
        new_tree(g, h_new, i_g + 1, i_h + 1)


def mcsp(g, h):
    global mcs_count
    g.append('_')
    h.append('_')
    mcsp_tree(g, h, 0, 0)
    print(mcs_count)


def new_search(g, h):
    global new_count
    # g_len = len(g)

    g_combinations = list(combinations(list(range(len(g))), 2))
    h_combinations = list(combinations(list(range(len(h))), 2))
    for (g1, g2) in g_combinations:
        g_new = g
        g_new = [g_new[g1]] + g_new[0:g1] + g_new[g1+1:len(g_new)]
        g_new = [g_new[g2]] + g_new[0:g2] + g_new[g2+1:len(g_new)]
        for (h1, h2) in h_combinations:
            h_new = h
            h_new = [h_new[h2]] + h_new[0:h2] + h_new[h2+1:len(h_new)]
            h_new = [h_new[h1]] + h_new[0:h1] + h_new[h1+1:len(h_new)]
            new_tree(g_new, h_new, 2, 2)
    print(new_count)


# g = [1, 2, 3, 4]
# h = ['a', 'b', 'c', 'd']
# print(list(combinations(list(range(len(g))), 2)))
g = list(range(6))
h = list(range(6))
new_search(g, h)
mcsp(g, h)

