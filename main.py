import sys
import utility, kmer, graph, quality_control

k = 20
path = ""
clean = True
test = False

"""
Attempt to get the proper arguments; kmer size, file path
"""
for arg in sys.argv:
    if arg.startswith("k="):
        k = int(arg[2:])
    if arg.startswith("f="):
        path = str(arg[2:])
    if arg.startswith("c="):
        clean = arg[2:]
    if arg.startswith("t="):
        test = arg[2:]

"""
Clean the data
"""
if clean == True:
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
g.save_graph()

"""
Save the Eulerian path (equivalent to the sequenced genome); print it
"""
seq = g.save_euler()
print("eulerian sequence length: %i" % len(seq))

"""
Play an audio cue so that we know the sequencing has completed; useful for long
runs when we are testing inefficient methods
"""
utility.play()
