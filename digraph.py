import sys, time
import kmer, sqldb
from Bio.Blast import NCBIWWW
from timer import timer
from graph import graph as gr

"""
This class is used to create an edge between two nodes in our graph
The label of the edge is the k-mer, or the paired k-mers after they've been
concatenated. It's stored in an SQL table with the edge index as a lookup number
"""
class edge:
    def __init__(self, db, index, label):
        self.db = db
        self.index = index
        self.multiplicity = 1

        if self.db.edges.has(label):
            self.index = self.db.edges.get_id(label)
        else:
            self.db.edges.add(index, label)

    def label(self):
        return self.db.edges.get_str(self.index)

    def __hash__(self):
        """ A custom hash so we can see if our edge is contained in a set """
        return hash(self.label)

    def __iter__(self):
        """ Allows us to iterate through the edge lists via name """
        return label

"""
This class is used to create a node on the graph, representing either the prefix
or suffix of the k-mer. For a set of paired k-mers, this is the prefix/suffix of
both k-mers, concatenated for the label. The prefix/suffix is stored in a lookup
table to cut memory usage. Stores connecting edge information
"""
class node:
    def __init__(self, db, index, label, paired=False):
        self.db = db
        self.edges = []
        self.index = index

        self.in_degree = 0
        self.out_degree = 0

        if self.db.nodes.has(label):
            self.index = self.db.nodes.get_id(label)
        else:
            self.db.nodes.add(index, label)

    """ return the label of the node from the database """
    def label(self):
        return self.db.nodes.get_str(self.index)

    def __eq__(self, other):
        return isinstance(other, self.__class__) and getattr(other,'label', None)

    def __hash__(self):
        """ A custom hash so we can see if our edge is contained in a set """
        return hash(self.db.get(index))

    def __str__(self):
        return "%s[%s]" % (self.index, self.db.get(index))

    def get_degree(self):
        return self.in_degree - self.out_degree
"""
This will create a directed graph (Debruijn graph) from a set of k-mers
Specify the k-length of the kmers, as well as the kmers previously created, and
whether or not this is a paired graph
"""
class graph:
    def __init__(self, k, paired=True, kmers=None):
        self.x = []
        self.test = True
        self.k = k
        self.db = sqldb.db("db/mydb")
        self.edges = dict()
        self.nodes = dict()
        self.indices = dict()
        self.contigs = dict()
        self.count = 0
        self.edge_index = 0

        self.paired = paired

        if kmers != None:
            self.add_kmers(kmers,paired)

            for i in kmers:
                e = edge(self.db, i)
                self.add_edge(e, paired)

    def add_kmers(self, km):
        for i,v in enumerate(self.labels_from_kmers(self.k, km)):
            if (self.nodes.get(v) == None):
                self.nodes[v] = node(self.db, self.count, v, self.paired)
                self.indices[self.count] = v
                self.count += 1

        for i in km:
            e = edge(self.db, self.edge_index, i)
            self.add_edge(e)
            self.edge_index += 1

    """ Create an edge, then add it to the graph; where the passed edge
    is a kmer string """
    def add_edge(self, edge):
        klen = self.k * 2
        hlen = int(klen/2)

        if not self.paired:
            lindex = edge.label()[:-1]
            rindex = edge.label()[1:]
        else:
            lindex = edge.label()[:hlen-2] + edge.label()[hlen:klen-1]
            rindex = edge.label()[1:hlen-1]+edge.label()[hlen+1:klen]

        """ Create the edges and update the degrees for the nodes """
        if edge.label() not in self.edges.keys():
            self.edges[edge.label()] = edge
            self.nodes[lindex].edges.append(edge)

            self.nodes[rindex].in_degree += 1
            self.nodes[lindex].out_degree += 1
        else:
            """ Currently unused; may be valid for determining low occurence
            rates for data cleaning """
            self.edges[edge.label()].multiplicity += 1



    """ Creates a sorted set of labels from the kmers to be used for node
    creation """
    def labels_from_kmers(self, k, kmers):
        if not self.paired:
            l = set([i[:k-1] for i in kmers] + [i[1:k] for i in kmers])
            return l
        else:
            klen = self.k*2 # cache length
            hlen = int(klen/2) # cache half length
            l = set([i[:hlen-2]+i[hlen:klen-1] for i in kmers]+[i[1:hlen-1]+i[hlen+1:klen] for i in kmers])

            return l


    """ Returns the Debruijn graph as a GraphViz source string """
    def export_graphviz(self):
        result = ''
        result += 'digraph {\n'
        result += '   graph [nodesep=2, size="300,300"];\n'
        for i in self.nodes:
            node = self.nodes[i]

            lenl = int(len(node.label())/2)
            l = str(node.index)+"["+node.label()[:lenl] + "_" + node.label()[lenl:]+"]"

            result += '    N%d [shape="box", style="rounded", label="%s"];\n' % (node.index, l)

        klen = self.k*2
        hlen = int(klen/2)

        for i in self.edges:
            if self.paired:
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
            if self.nodes[n].in_degree == 0:
                return self.nodes[n]

    """ When called by get_chained_contigs, this will return the path travelled until
    it reaches a branching point with out degree < 1 """
    def get_contig_region(self,start = None, next = None):
        node = start

        """"
        If we specified a node, set path = [start] and the node to check as next
        """
        if next == None:
            path = []
        else:
            node = next
            path = [start]

        """
        If path[start] != None, that means the edge (and trailing contigs)
        have already been dealt with
        """
        if len(path) > 0 and self.contigs[path[0].label] != None:
            return False

        klen = self.k*2
        hlen = int(klen/2)

        while node.out_degree <= 1:
            path.append(node)
            if (len(node.edges) > 0):
                node = self.nodes[node.edges[0].label()[1:hlen-1]+node.edges[0].label()[hlen+1:klen]]
            else:
                self.contigs[path[0].label()] = path
                return path[-1]

    """ Returns a list of all the contigs found starting from the given node;
    if no node is passed, the get_start() method is called to find a proper
    start """
    def get_chained_contigs(self, node=None):
        # Start at the beginning if we haven't started yet
        if node == None:
            node = self.get_start()

        node = self.get_contig_region(node)
        contigs = []

        klen = self.k*2
        hlen = int(klen/2)

        # If we've hit the end of the graph, start collapsing, else, continue
        while node != False and node.out_degree > 0:
            for edge in node.edges:
                if not paired:
                    lindex = edge.label()[:-1]
                    rindex = edge.label()[1:]
                else:
                    lindex = edge.label()[:hlen-2] + edge.label()[hlen:klen-1]
                    rindex = edge.label()[1:hlen-1]+edge.label()[hlen+1:klen]

                node = self.get_contig_region(self.nodes[lindex],self.nodes[rindex])


    """ Turns the contigs in to their textual representations (with the kmers
    concatenated as a superstring) and saves them to a file.
    Default save location: 'Data/contigs.dat'
    """
    def save_contigs(self, path="Data/contigs.dat"):
        self.get_chained_contigs()

        start = True
        hlen = int(self.k/2)
        contigs = self.contigs

        with open(path, 'w') as file:
            for contig in contigs:
                for node in contigs[contig]:
                    if start:
                        file.write(node.label()[:hlen-1])
                        start = False
                    file.write(node.label()[:hlen][-1])


        print("[DONE]   ->  Contiguous regions saved; %i found"%len(contigs))

    """ Gets the GraphViz source from the export function, then saves it in a
    file to be used later by a program like gvedit, or the gv.save() function """
    def save_graph(self):
        file = open("output.gv","w")
        file.write(self.export_graphviz())
        file.close()
        print("[DONE]   ->  Saving output as GraphViz compatible format")

if __name__ == "__main__":
    k = 31

    start = time.time()

    g = graph(k,True)
    kmer.graph_from_sequences(g, k, "Data/lseq.dat", "Data/rseq.dat")

    g.save_contigs()
    g.save_graph()

    """
    print("[TASK]   ->  Performing basic local alignment search for resulting nucleotide sequence")
    sys.stdout.flush()

    nt_string = open("Data/contigs.dat").read()
    result_handle = NCBIWWW.qblast("blastn", "nt", nt_string)

    sys.stdout.flush()
    with open("blast_results.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
        result_handle.close()


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
    """
    print("Elapsed: %f"%(time.time() - start))
