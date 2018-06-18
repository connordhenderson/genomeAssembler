import sys, ast, random

def random_sequence(read_length, insert, read_amt, path):
    count = 0
    seq = ''
    llist = []
    rlist = []

    with open(path) as file:
        try:
            while 1:
                line = next(file)
                count += 1
                seq += line
        except StopIteration:
            seq_length = len(seq)
            print("scanned %s; %i lines"%(path,count))
            print("seq length: %i"%(seq_length))

    while read_amt > 0:
        index = random.randint(0,seq_length)
        rindex = index + insert + read_length
        if (seq_length > index+2*read_length+insert):
            l = seq[index:index+read_length]
            r = seq[rindex:rindex+read_length]
            llist.append(l)
            rlist.append(reverse_compliment(r))
            read_amt -= 1
    with open("Data/lseq.dat", 'w') as lfile, open("Data/rseq.dat", 'w') as rfile:
        for x in llist:
            lfile.write("%s\n"%x)
        for line in rlist:
            rfile.write("%s\n"%line)
    print("[DONE]   ->  paired sequences created & saved")
    return seq

for argv in sys.argv:
    if argv.startswith("f="):
        path = argv[2:]
    if argv.startswith("d="):
        d = ast.literal_eval(argv[2:])

def reverse_compliment(seq):
    seq = list(seq)
    compliment = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    for i in range(len(seq)):
        seq[i] = compliment[seq[i]]
    return "".join(reversed(seq))

random_sequence(d[0],d[1],d[2],path)
