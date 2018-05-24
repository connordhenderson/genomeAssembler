import os, random

def create_kmers(k, inpath="Data/sequences.dat"):
    old_kmers = set()
    with open(inpath) as handle:
        try:
            while 1:
                line = next(handle)
                if (line != "\n"):
                    new_kmers = kmer_from_string(line, k)
                    for kmer in new_kmers:
                        old_kmers.add(kmer)

        except StopIteration:
            print("[DONE]   ->  kmers created")
    return sorted(list(old_kmers))



def kmer_from_string(string, k):
    k = int(k)
    string = string.strip()
    kmers = sorted([string[i:i+k] for i in range(len(string) - k+1)])
    """
    # Originally used for allowing duplications; no longer necessary since we
    # are relying on creating contigs first
    for i in range(1,len(kmers)):
        if kmers[i] == kmers[i-1][0:k]:
            t = kmers[i-1].find('_')
            if (t >= 0):
                # TODO: Address the route of the ValueError
                n = int(kmers[i-1][t+1:]) + 1
                kmers[i] = kmers[i] + "_" + str(n)
            else:
                kmers[i-1] = kmers[i-1] + "_1"
                kmers[i] = kmers[i] + "_2"
    """
    return kmers

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
