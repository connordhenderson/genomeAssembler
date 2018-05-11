import sys
import time
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

fpath = sys.argv[1]

# Get PHRED score
pqual = 20
if (len(sys.argv) >= 3):
    pqual = sys.argv[2]


seqrec = []

# Filter our data by a given PHRED score
start = time.time()

count = 0
h = []
with open(fpath) as handle:
    for title, seq, qual in FastqGeneralIterator(handle):
        count += 1
        if (min([ord(i) for i in qual]) <= 20):
            h.append(seq)

print ("records: %i   -> elapsed: %f" % (count, (time.time() - start)))

"""
records = (rec for rec in SeqIO.parse(open(fpath),"fastq")
            if min(rec.letter_annotations["phred_quality"]) <= pqual)

handle = open("Data/trimmed1.fastq", "w")
count = SeqIO.write(records, handle, "fastq")
"""
