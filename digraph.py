import sys
import kmer
from Bio.Blast import NCBIWWW

"""
Creates an edge between two nodes in our graph
The label of the edge is the k-mer, or the paired k-mers after they've been
concatenated
"""
class edge:
    def __init__(self, label):
        self.label = label
        self.multiplicity = 1

    def __hash__(self):
        """ A custom hash so we can see if our edge is contained in a set """
        return hash(self.label)

    def __iter__(self):
        """ Allows us to iterate through the edge lists via name """
        return label

"""
Creates a node on the graph, representing either the prefix or suffix of the
k-mer. For a set of paired k-mers, this is the prefix/suffix of both k-mers,
concatenated for the label. Stores connecting edge information
"""
class node:
    def __init__(self, graph, index, label, paired=False):
        self.graph = graph
        self.label = label
        self.edges = []
        self.index = index

        self.in_degree = 0
        self.out_degree = 0

    def __eq__(self, other):
        return isinstance(other, self.__class__) and getattr(other,'label', None)

    def __hash__(self):
        """ A custom hash so we can see if our edge is contained in a set """
        return hash(self.label)

    def __str__(self):
        return "%i[%s]" % (self.index, self.label)

    def get_degree(self):
        return self.in_degree - self.out_degree
"""
This will create a directed graph (Debruijn graph) from a set of k-mers
Specify the k-length of the kmers, as well as the kmers previously created, and
whether or not this is a paired graph
TODO: Verify unpaired graph creation is unaffected by implementation of pairs
"""
class graph:
    def __init__(self, k, kmers, paired=False):
        self.k = k
        self.edges = dict()
        self.nodes = dict()
        self.indices = dict()
        self.contigs = dict()

        for i,v in enumerate(self.labels_from_kmers(k,kmers, paired)):
            self.nodes[v] = node(self,i, v, paired)
            self.indices[i] = v

        print("[DONE]   ->  Edges created")

        for i in kmers:
            e = edge(i)
            self.add_edge(e, paired)

    """ Create an edge, then add it to the graph; where the passed edge
    is a kmer string """
    def add_edge(self, edge, paired=False):
        klen = len(kmers[0])
        hlen = int(klen/2)

        if not paired:
            lindex = edge.label[:-1]
            rindex = edge.label[1:]
        else:
            lindex = edge.label[:hlen-2] + edge.label[hlen:klen-1]
            rindex = edge.label[1:hlen-1]+edge.label[hlen+1:klen]

        """ Create the edges and update the degrees for the nodes """
        if edge.label not in self.edges.keys():
            self.edges[edge.label] = edge
            self.nodes[lindex].edges.append(edge)
        else:
            """ Currently unused; may be valid for determining low occurence
            rates for data cleaning """
            self.edges[edge.label].multiplicity += 1

        self.nodes[rindex].in_degree += 1
        self.nodes[lindex].out_degree += 1

    """ Creates a sorted set of labels from the kmers to be used for node
    creation """
    def labels_from_kmers(self, k, kmers, paired=False):
        if not paired:
            l = sorted(set([i[:k-1] for i in kmers] + [i[1:k] for i in kmers]))
            return l
        else:
            klen = len(kmers[0]) # cache length
            hlen = int(klen/2) # cache half length
            l = sorted(set([i[:hlen-2]+i[hlen:klen-1] for i in kmers]+[i[1:hlen-1]+i[hlen+1:klen] for i in kmers]))

            return l


    """ Returns the Debruijn graph as a GraphViz source string """
    def export_graphviz(self, paired=False):
        result = ''
        result += 'digraph {\n'
        result += '   graph [nodesep=2, size="300,300"];\n'
        for i in self.nodes:
            node = self.nodes[i]

            lenl = int(len(node.label)/2)
            l = str(node.index)+"["+node.label[:lenl] + "_" + node.label[lenl:]+"]"

            result += '    N%d [shape="box", style="rounded", label="%s"];\n' % (node.index, l)

        klen = self.k*2
        hlen = int(klen/2)

        for i in self.edges:
            if paired:
                src = self.nodes[i[:hlen-2] + i[hlen:klen-1]].index
                dst = self.nodes[i[1:hlen-1] + i[hlen+1:klen]].index
            else:
                src = self.nodes[i[:-1]].index
                dst = self.nodes[i[1:]].index
            result += '    N%d -> N%d' % (src, dst)
            if (len(i) > 0):
                result += ' [label="%s", penwidth=1.0]' % i
            result += ';\n'
        result += '    overlap=false;\n'
        result += '}\n'
        return result

    """ Returns the first node found with an in degree of 0; theoretically
    this should be the starting point """
    # TODO: Remove nodes/edges in contig from list, find more contigs
    def get_start(self):
        for n in self.nodes:
            if self.nodes[n].get_degree() == -1 and self.nodes[n].in_degree == 0:
                return self.nodes[n]

    """ Returns a list of all the contigs found starting from the given node;
    if no node is passed, the get_start() method is called to find a proper
    start """
    def get_all_contigs(self, paired=False, node=None):
        # Start at the beginning if we haven't started yet
        if node == None:
            node = self.get_start()

        node = self.get_contig(node)
        contigs = []

        klen = self.k*2
        hlen = int(klen/2)

        # If we've hit the end of the graph, start collapsing, else, continue
        while node != False and node.out_degree > 0:
            for edge in node.edges:
                if not paired:
                    lindex = edge.label[:-1]
                    rindex = edge.label[1:]
                else:
                    lindex = edge.label[:hlen-2] + edge.label[hlen:klen-1]
                    rindex = edge.label[1:hlen-1]+edge.label[hlen+1:klen]

                node = self.get_contig(self.nodes[lindex],self.nodes[rindex])


    """ When called by get_all_contigs, this will return the path travelled until
    it reaches a branching point with out degree < 1 """
    def get_contig(self,start = None, next = None, paired=False):
        node = start
        if next == None:
            path = []
        else:
            node = next
            path = [start]

        klen = self.k*2
        hlen = int(klen/2)

        while node.out_degree <= 1:
            path.append(node)
            if (len(node.edges) > 0):
                node = self.nodes[node.edges[0].label[1:hlen-1]+node.edges[0].label[hlen+1:klen]]
            else:
                self.contigs[path[0].label] = path
                return path[-1]

    """ Turns the contigs in to their textual representations (with the kmers
    concatenated as a superstring) and saves them to a file.
    Default save location: 'Data/contigs.dat'
    """
    def save_contigs(self, paired=False, path="Data/contigs.dat"):
        self.get_all_contigs(paired)

        start = True
        hlen = int(self.k/2)
        contigs = self.contigs

        with open(path, 'w') as file:
            for contig in contigs:
                for node in contigs[contig]:
                    if start:
                        file.write(node.label[:hlen-1])
                        start = False
                    file.write(node.label[:hlen][-1])


        print("[DONE]   ->  Contiguous regions saved; %i found"%len(contigs))

    """ Gets the GraphViz source from the export function, then saves it in a
    file to be used later by a program like gvedit, or the gv.save() function """
    def save_graph(self, paired=False):
        file = open("output.gv","w")
        file.write(self.export_graphviz(paired))
        file.close()
        print("[DONE]   ->  Saving output as GraphViz compatible format")


k = 31
kmers = kmer.create_kmers(k,"Data/lseq.dat","Data/rseq.dat")

g = graph(k,kmers,True)

g.save_contigs(True)
g.save_graph(True)

"""
print("[TASK]   ->  Performing basic local alignment search for resulting nucleotide sequence")
sys.stdout.flush()

nt_string = open("Data/contigs.dat").read()
result_handle = NCBIWWW.qblast("blastn", "nt", nt_string)

sys.stdout.flush()
with open("blast_results.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
    result_handle.close()
"""

print("[DONE]   ->  NCBI nucleotide blast top match:")
from xml.dom import minidom

def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
    return ''.join(rc)

xmldoc = minidom.parse("blast_results.xml")
hitlist = xmldoc.getElementsByTagName('Hit')

count = 0
for node in hitlist:
    if count < 0:
        break
    count -= 1

    hit_ids = node.getElementsByTagName('Hit_id')
    for id in hit_ids:
        print(id.childNodes[0].nodeValue)

    Hit_def = node.getElementsByTagName('Hit_def')
    for id in Hit_def:
        print(id.childNodes[0].nodeValue)
