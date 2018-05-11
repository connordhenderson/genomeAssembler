import sys, itertools

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def find_adapter(path, n, length):
    # Looks through the first N records to try and find any adapters
    records = []
    _n = n

    print("getting %i records of length %i" % (n, length))

    with open(path) as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            if (int(n) <= 0):
                break
            else:
                records.append(seq[:length])
                n -= 1
    #for rec in records:
    #    print(rec)
    print("records: %i" % len(records))

    matches = 0
    print (range(_n))
    for i in range(_n):
        for j in range(_n):
            if i != j and j <_n-matches and i<_n-matches:
                if records[i] == records[j]:
                    print (records[i])
                    records.pop(i)
                    matches += 1

print("starting")
find_adapter(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
