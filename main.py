import sys
import time
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGIter




if (len(sys.argv) < 3):
    print("Wrong args")
    sys.exit()
fpath = sys.argv[1]
pqual = sys.argv[2]

# we want to trim the data of anything that falls below our PHRED quality
qc_data = (rec for rec in SeqIO.parse(fpath,"fastq")
            if min(rec.letter_annotations["phred_quality"]))

for rec in qc_data:
    print(rec.seq)
