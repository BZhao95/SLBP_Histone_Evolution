import nuad, nupack, scadnano
import nuad.vienna_nupack
import stem_loops as sl
# from stem_loops import get_energy, get_stem1, get_loop, get_stem2

'''
Library for computing distances between stem loops

Example functions, each take stem loops as input:
 - compute edit distances between two stem loops
 - compute an energy distance between two stem loop;  use different energy models (e.g., NUPACK, ViennaRNA) to compute energy distances
'''

def DG_distance(sl1:str, sl2:str, temperature:float=27, method:str='nupack mfe') -> float:
    return sl.get_energy(sl1.replace(" ",""), temperature, method) - get_energy(sl2.replace(" ",""), temperature, method)

def stem_edit_distance(sl1:str, sl2:str) -> int:
    '''number of base pair mutations needed to convert stem of sl1 into stem of sl2
    sl1, sl2 are assumed to be in domain format with domains separated by spaces'''
    s1 = sl.get_stem1(sl1)
    s2 = sl.get_stem1(sl2)
    return sum([1 for i in range(len(s1)) if s1[i] != s2[i]])

def loop_edit_distance(sl1:str, sl2:str) -> int:
    '''number of base pair mutations needed to convert loop of sl1 toloop of sl2'''
    l1 = sl.get_loop(sl1)
    l2 = sl.get_loop(sl2)
    return sum([1 for i in range(len(l1)) if l1[i] != l2[i]])

def stem_loop_edit_distance(sl1:str, sl2:str) -> int:
    '''number of base pair edits needed to convert stem and loop of sl1 into sl2'''
    return stem_edit_distance(sl1, sl2) + loop_edit_distance(sl1, sl2)

# def build_tree_from_stem_loops(stem_loops:List[str], distance_func:Callable[[str,str],float]) -> nx.Graph:
#     '''build a tree from stem loops using a distance function to compute edge weights'''
#     G = nx.Graph()
#     for i in range(len(stem_loops)):
#         for j in range(i+1, len(stem_loops)):
#             G.add_edge(i, j, weight=distance_func(stem_loops[i], stem_loops[j]))
#     return G


import networkx as nx
from typing import List, Callable

def build_graph_from_stem_loops(stem_loops:list[str], distance_func:Callable[[str,str],float]) -> nx.Graph:
    '''build a graph from stem loops using a distance function to compute edge weights'''
    G = nx.Graph()
    for i in range(len(stem_loops)):
        for j in range(i+1, len(stem_loops)):
            if distance_func(stem_loops[i], stem_loops[j]) == 1:
                G.add_edge(i, j, weight=distance_func(stem_loops[i], stem_loops[j]), )
    return G