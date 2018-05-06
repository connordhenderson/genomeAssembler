import sys
import itertools
import time

# Given a list of strings 'str' and a size of a k-mer 'k', this will construct the set of nodes/edges for a De Bruijn graph
def create_graph(str,k):
    edges = []
    nodes = set()
    for i in range(len(str) - k + 1):
        # Creating left/right k-mers
        edges.append([str[i:i+k-1], str[i+1:i+k]])
        nodes.add(str[i:i+k-1])
        nodes.add(str[i+1:i+k])
    return nodes, edges

# Determines how much of the end of 'a' overlaps with the beginning of 'b'
def overlap(a, b, min_length = 3):
    start = 0
    while 1:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

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

def get_eulerian(path):
    return list(itertools.permutations(path,2))

seq = "ACGCGTCG"
nodes, edges = create_graph(seq, 4)
print ("seq:   %s"%seq)
print ("nodes: %s"%nodes)
print ("edges: %s"%edges)


start = time.time()

st = "to_everything_turn_turn_turn_there_is_season"
G = create_graph(st,4)
path = get_eulerian(G[1])

gen_dict = dict()
for rec in path:
    print(rec[0])
    gen_dict[rec[0]] = rec[1]

print(gen_dict["to_"])

print("elapsed: %f"% (time.time() - start))
