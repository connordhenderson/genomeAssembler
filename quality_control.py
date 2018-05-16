from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

"""
Trims the reads based on their Illumina score converted to a Phred quality
Illumina quality scores have an offset of 33 compared to their Phred score
"""
def trim_quality(reads, pqual):
    for read in reads:
        if record.letter_annotations["phred_quality"] >= pqual:
            yield read

"""
Given a primer, create a generator which will return the sequences with the
primer trimmed from the beginning
"""
def trim_primer(reads, primer):
    len_primer = len(primer)
    for read in reads:
        if record.seq.startswith(primer):
            yield record[len_primer:]
        else:
            yield record

"""
Given a list of adapters, create a generator which will only return the
sequence that comes after the adapter
"""
def trim_adapters(reads, adapter):
    len_adapter = len(adapter)
    for read in reads:
        index = read.seq.find(adapter)
        if index == -1:
            yield read
        else:
            yield read[index+len_adapter:]

"""
Performs all of the various data trimming specified
"""
def trim_data(reads, pqual=20, primer=None, adapter=None):
    pqual += 33
    reads = trim_quality(reads,pqual)
    if primer != None:
        reads = trim_primer(reads,primer)
    if adapter != None:
        reads = trim_adapter(reads,adapter)
    return reads

"""
Initiates the process to trim the data, then write it to a file
"""
def clean_data(path):
    reads = SeqIO.parse(path,"fastq")
    reads = trim_data(reads)

    print("[DONE]   ->  data cleaning generators created... saving sequences")
    seqlist = []
    with open(path) as handle:
        for title,seq,qual in FastqGeneralIterator(handle):
            seqlist.append(seq)

    path = "Data/sequences.dat"

    file = open(path, 'w')
    for item in seqlist:
        file.write("%s\n" % item)
