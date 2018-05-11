import sys, time
import utility as util
import graph, gnode
import itertools
import time

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from collections import defaultdict

# Finds the strings with the highest overlap values, and returns those values
def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_overlap_len = 0
    for a,b in itertools.permutations(reads, 2):
        overlap_len = overlap(a,b,min_length = k)
        if overlap_len > best_overlap_len:
            reada, readb = a,b
            best_overlap_len = overlap_len
    return reada, readb, best_overlap_len

# Uses the maximal overlaps to always combine the highest overlap values in order; not exhaustive but not always right.
# This means we can get stuck in a cycle and need to find a way to break out of the cycle
def greedy_shortest_common_superstring(reads, k):
    reada, readb, overlap_len = pick_maximal_overlap(reads, k)
    while overlap_len > 0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada + readb[overlap_len:])
        reada, readb, overlap_len = pick_maximal_overlap(reads, k)

    return ''.join(reads)



# This allows us to have duplicate keys since they're all assigned via ID
def number_vertices(edges):
    v = dict()

    curr_num = 0;
    for e in edges:
        if (e[0] not in v):
            v[e[0]] = curr_num
            curr_num += 1

    # We actually want to reverse the order of the dict for looking up, and
    # eliminating duplicates
    r = dict()
    for n in list(v):
        r[v[n]] = n

    return v, r

def create_numbered_dict(edges, lookup):

    d = dict()
    for e in edges:
        if (lookup[e[0]] in d):
            arr = []
            arr.append(d[lookup[e[0]]])
            arr.append(lookup[e[1]])
            d[lookup[e[0]]] = arr
        else:
            d[lookup[e[0]]] = lookup[e[1]]
    return d

sequence1 = "this_is_my_test_sentence_for_graphing"
sequence2 = "this_branches_from_previous_sequence!"
nodes, edges = graph.create_graph(None, None, sequence1, 6)
nodes, edges = graph.create_graph(nodes, edges, sequence2, 6)



#edges = test(seqrec,3)

#print(len(edges), len(seqrec))
