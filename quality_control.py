from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

"""
Trims the reads based on their Illumina score converted to a Phred quality
Illumina quality scores have an offset of 33 compared to their Phred score
"""
def trim_quality(reads, pqual):
    for read in reads:
        if min(record.letter_annotations["phred_quality"]) >= pqual:
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
def trim_data(reads, pqual=640, primer=None, adapter=None):
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
def clean_data(path, save):
    reads = SeqIO.parse(path,"fastq")
    reads = trim_data(reads)

    seqlist = []

    good_reads = (rec for rec in \
              SeqIO.parse(path, "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 20)
    for rec in good_reads:
        seqlist.append(rec.seq)



    file = open(save, 'w')
    print("saving to '%s'"%save)
    for item in seqlist:
        file.write("%s\n" % item)

"""
Initiates the process to trim the data, then write it to a file
"""
def clean_paired_data(lpath, rpath, lsave, rsave):
    pqual = 40

    lseqlist = []
    rseqlist = []
    count = 0
    with open(lpath) as lfile, open(rpath) as rfile:
        liter = FastqGeneralIterator(lfile)
        riter = FastqGeneralIterator(rfile)

        while 1:
            try:
                lseq, lqual = next(liter)[1:]
                rseq, rqual = next(liter)[1:]

                if ord(min(lqual)) >= pqual and ord(min(rqual)) >= pqual:
                    lseqlist.append(lseq)
                    rseqlist.append(rseq)

            except StopIteration:
                print("EOF")
                break

    with open(lsave,'w') as lfile, open(rsave, 'w') as rfile:
        print("saving to '%s' and '%s'"%(lsave,rsave))
        for litem in lseqlist:
            lfile.write("%s\n" % litem)
        for ritem in rseqlist:
            rfile.write("%s\n" % ritem)
