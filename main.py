import sys, time
import utility, quality_control

def start(test=False):
    start = time.time()

    k = 20
    lpath = ""
    rpath = ""
    clean = True
    test = False

    """
    Attempt to get the proper arguments; kmer size, file path
    """
    for arg in sys.argv:
        if arg.startswith("k="):
            k = int(arg[2:])
        if arg.startswith("l="):
            lpath = str(arg[2:])
        if arg.startswith("r="):
            rpath = str(arg[2:])
        if arg.startswith("t="):
            test = arg[2:]

    """
    Clean the data
    """
    if test == False:
        quality_control.clean_paired_data(lpath, rpath, "Data/lseq.dat", "Data/rseq.dat")

        print("[DONE]   ->  data trimmed")
    else:
        print("[SKIP]   ->  data cleaning skipped")


if sys.argv[0] == 'main.py':
    start()
