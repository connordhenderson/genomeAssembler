import sys, time
import utility, kmer, graph, quality_control
def start(test):
    start = time.time()

    k = 20
    path = ""
    clean = True

    """
    Attempt to get the proper arguments; kmer size, file path
    """
    for arg in sys.argv:
        if arg.startswith("k="):
            k = int(arg[2:])
        if arg.startswith("f="):
            path = str(arg[2:])
        if arg.startswith("t="):
            test = arg[2:]

    """
    Clean the data
    """
    if test == False:
        quality_control.clean_data(path)
        print("[DONE]   ->  data trimmed")
    else:
        print("[SKIP]   ->  data cleaning skipped")

    """
    Create the k-mers
    """
    kmer.clear_kmers("Data/kmers.dat")
    kmers = kmer.create_kmers(k)

    print("# of kmers: %i" % len(kmers))
    sys.stdout.flush()      # Flush stdout so we can identify when we get stuck in
                            # finding the Eulerian tour
    """
    Create the graph from the k-mers
    """
    g = graph.create_graph(kmers, k)
    sys.stdout.flush()

    if len(kmers) < 5000:
        g.save_graph()
        sys.stdout.flush()
    else :
        print("[SKIP]   ->  The graph is too large to display; skipping output.gv")
        sys.stdout.flush()


    """
    Save the Eulerian path (equivalent to the sequenced genome); print it
    """
    print("[TASK]   ->  Finding Eulerian Path")
    sys.stdout.flush()
    seq = g.save_euler()
    if (seq != False):
        print("eulerian sequence length: %i" % len(seq))

    sys.stdout.flush()
    print("time elapsed: %f" % (time.time() - start))

if sys.argv[0] == 'main.py':
    start()
