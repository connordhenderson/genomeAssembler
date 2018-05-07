import sys
import itertools
import time

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from collections import defaultdict
from drawDebruijn import draw as drawDeb

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

def test(reads, k):
    read_dict = dict()
    list = []


    index = 0
    for a in reads:
        current = None
        paired = None

        # We want to remove this element so we aren't comparing it against itself
        reads.pop(index)

        for b in reads:
            overlap_amt = overlap(a,b, min_length = 1)
            if current == None or current < overlap_amt:
                current = overlap_amt
                paired = b
            # We now have the highest overlap for a read
            read_dict[a] = b

        reads.insert(index,a)

        # found the best match
        list.append([a,paired])
        index += 1
        current = None

    return list


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

# We're going to attempt to implement Hierholzer's algorithm
def get_eulerian(edges):
    end = False

    d = create_dict(edges)

    # we need to create a dictionary of our edges
    for e in edges:
        d[e[0]] = e[1]

    v = list(d)[0]
    stack = [v]
    circuit = []

    current = 0
    next = 0

    # count the amount of edges on each vertex
    e_count = dict()
    for vx in range(len(list(d))):
        e_count[vx] = len(d[vx])

    while len(c_path) >= 0:
        # if we have an exit edge, carry on
        if len(d[vx]) > 0:
            stack.append(current)

            next = d[current][0]
            # remove that edge
            arr = d[current]
            arr = arr[1:]
            d[current] = arr

            current = next
        else:
            circuit.append(current)
            current = stack[-1]
            stack.pop()
    return circuit


    # Choose a starting vertex v, follow edges until we hit a closed loop at v


# We need to be able to create a dictionary that will store information
def create_dict(edges):
    d = dict()
    for e in edges:
        if (e[0] in d):
            arr = []
            arr.append(d[e[0]])
            arr.append(e[1])
            d[e[0]] = arr
        else:
            d[e[0]] = [e[1]]
    return d

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

def print_circuit(circuit):
    return "d"


good_data = (rec for rec in SeqIO.parse(sys.argv[1],"fastq")
            if min(rec.letter_annotations["phred_quality"]) >= 20)

seqrec = []
for x in good_data:
    seqrec.append(str(x.seq))

edges = test(seqrec, 3)
print(edges)

drawDeb(edges)
#edges = test(seqrec,3)

#print(len(edges), len(seqrec))
