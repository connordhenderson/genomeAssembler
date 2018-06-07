import os, random, sys
from digraph import graph

class kmer:
    def __init__(self, l, r=""):
        self.seq = l + r

    def startswith(self, kmer):
        return self.l.startswith(kmer.l)

    def __getitem__(self, key):
        return self.seq[key]

    def __iter__(self):
        return self.seq

    def __eq__(self, s):
        return self.seq == s

    def __lt__(self, other):
        return self.seq < other.seq

    # We want to hash by the sequence for set membership, not comparing ID
    def __hash__(self):
        return hash(self.seq)

    def __str__(self):
        return self.seq

def create_kmers(k, lpath="Data/sequences.dat", rpath = None):
    old_kmers = []
    id = 0
    print("DEPREC")
    """
    if rpath == None:
        with open(lpath) as handle:
            try:
                while 1:
                    line = next(handle)
                    if (line != "\n"):
                        new_kmers = kmer_from_string(line, id, k)
                        for km in new_kmers:
                            if km not in old_kmers:
                                old_kmers[km] = kmer(id, km)
                            else:
                                old_kmers[km].id.append(id)
                            id += 1
            except StopIteration:
                print("[DONE]   ->  kmers created")
    else:
    """
    with open(lpath) as lhandle, open(rpath) as rhandle:
        try:
            while 1:
                lline = next(lhandle)
                rline = next(rhandle)

                if (lline != "\n" and rline != "\n"):
                    new_kmers = paired_kmer_from_string(lline, rline, id, k)
                    for km in new_kmers:
                        if km not in old_kmers:
                            old_kmers.append(km)
        except StopIteration:
            print("[DONE]   ->  K-mers created")
    return sorted(list(old_kmers))

def graph_from_sequences(graph, k, lpath="Data/sequences.dat", rpath = None):
    old_kmers = dict()
    """
    if rpath == None:
        with open(lpath) as handle:
            try:
                while 1:
                    line = next(handle)
                    if (line != "\n"):
                        new_kmers = kmer_from_string(line, id, k)
                        graph.add_kmers(new_kmers)
            except StopIteration:
                print("[DONE]   ->  kmers created")
    else:
    """
    with open(lpath) as lhandle, open(rpath) as rhandle:
        try:
            print("[TASK]   ->  Creating k-mers")
            sys.stdout.flush()

            while 1:
                lline = next(lhandle)
                rline = reverse_compliment(next(rhandle).rstrip('\n'))

                if (lline != "\n" and rline != "\n"):
                    new_kmers = paired_kmer_from_string(lline, rline, k)
                    graph.add_kmers(new_kmers)
                    new_kmers = None
        except StopIteration:
            print("[DONE]   ->  Paired k-mers created")

def reverse_compliment(seq):
    seq = list(seq)
    compliment = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    for i in range(len(seq)):
        seq[i] = compliment[seq[i]]
    return "".join(reversed(seq))

def paired_kmer_from_string(lstring, rstring, k):
    k = int(k)
    lstring = lstring.rstrip('\n')
    rstring = rstring.rstrip('\n')

    return [lstring[i:i+k]+rstring[i:i+k] for i in range(len(lstring) - k+1)]

def clear_kmers(path="Data/kmers.dat"):
    file = open(path, 'w')
    file.close()

def list_to_file(path, kmers):
    contents = []

    if not os.path.exists(path):
        file = open(path, 'w')
        file.close()
    else:
        file = open(path, 'r')
        contents = file.readlines()
        file.close()

    with open(path, 'w') as file:
        p_kmer = ''
        for kmer in sorted(contents + kmers):
            if (kmer.strip() != p_kmer.strip()):
                if kmer != "\n":
                    p_kmer = kmer
                    file.write(kmer.strip() + "\n")
        file.close()
