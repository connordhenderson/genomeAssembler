import os, sys, random, ast, time
import kmer, quality_control, graph, utility, main, gv

def clean_sample(path):
    lines = []
    with open(path) as file:
        count = 0
        try:
            while 1:
                count += 1
                line = next(file).strip()
                if count == 2:
                    lines.append(line)
                if count == 4:
                    count = 0
        except StopIteration:
            print("done")
        file.close()

    with open("Data/sequences.dat", 'w') as file:
        for line in lines:
            file.write("%s\n"%line)

def lorem_ipsum():
    seq = ''
    with open("Data/loremipsum.dat") as file:
        try:
            while 1:
                line = next(file)
                line = "".join(line.splitlines())
                seq += "".join(line)
        except StopIteration:
            print("done parsing lorem ipsum")
            file.close()

    with open("Data/sequences.dat", 'w') as file:
        file.write("%s\n"%seq)
        file.close()


def random_sequence(seq_length, read_length, read_amt, path=None):
    seq_list = []

    chars = ["A","T","G","C"]
    seq = ''

    if path == None:
        for i in range(seq_length):
            seq += random.choice(chars)
    else:
        count = 0
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
        if (seq_length > index+read_length):
            s = seq[index:index+read_length]
            seq_list.append(s)
            read_amt -= 1
    with open("Data/sequences.dat", 'w') as file:
        for line in seq_list:
            file.write("%s\n"%line)
    print("[DONE]   ->  sequence created & saved")
    return seq

def create_paired_data(read_length, read_amt, gap_length, path):
    total_length = 2*read_length + gap_length
    seq = ''

    l_read = []
    r_read = []

    # Store the sequence in memory
    with open(path) as handle:
        try:
            while 1:
                seq += next(file)
        except StopIteration:
            pass

    seq_len = len(seq)

    while read_amt > 0:
        index = random.randint(0,len(seq))

        if seq_len > index + read_length:
            l_read.append(seq[index:index+read_length])
        if seq_len > index + total_length:
            new_index = index + gap_length + read_length
            r_read.append(seq[new_index:new_index + read_length])
        else:
            r_read.append(seq[-read_length:])

    with open("Data/left.dat", 'w') as file:
        for line in l_read:
            file.write("%s\n"%line)

    with open("Data/right.dat", 'w') as file:
        for line in r_read:
            new_line = ''
            for letter in reversed(line):
                # Get the reverse compliment (to simulate paired read)
                if letter == 'A':
                    new_line += 'T'
                if letter == 'T':
                    new_line += 'A'
                if letter == 'C':
                    new_line += 'G'
                if letter == 'G':
                    new_line += 'C'

            file.write("%s\n"%new_line)


print("Starting!")
sys.stdout.flush()

k = 10
path = None
test_data = None
for arg in sys.argv:
    if arg.startswith("k="):
        k = int(arg[2:])
    if arg.startswith("f="):
        path = str(arg[2:])
    if arg.startswith("d="):
        test_data = ast.literal_eval(arg[2:])
    if arg.startswith("l="):
        if arg[2:] == "create":
            seq = lorem_ipsum()
            sys.exit()
        if arg[2:] == "full":
            seq = lorem_ipsum()
            seq = random_sequence(test_data[0], test_data[1], test_data[2], path)
            main.start(True)
            #utility.play()
            sys.exit()
        if arg[2:] == "gbff":
            start = time.time()
            seq = random_sequence(test_data[0], test_data[1], test_data[2], path)
            main.start(True)
            #utility.play()
            print("[DONE]   ->  Elapsed time: %ss" % str(time.time() - start))

            sys.exit()

if (test_data != None) and len(test_data) > 2:
    print ("[TEST]  ->  random_sequence(",test_data,")")
    seq = random_sequence(test_data[0], test_data[1], test_data[2], path)
