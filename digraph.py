import sys
import kmer

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


class node:
    def __init__(self, graph, label, index):
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
        return "%i: %s" % (self.index, self.label)

    def get_degree(self):
        return self.in_degree - self.out_degree

class graph:
    def __init__(self, k, kmers):
        self.k = k
        self.edges = dict()
        self.nodes = dict()
        self.indices = dict()
        self.contigs = dict()

        for i,v in enumerate(self.labels_from_kmers(k,kmers)):
            self.nodes[v] = node(self,v,i)
            self.indices[i] = v

        """ Create the edges of the graph from our kmers """

        """
        l = []
        for i in range(len(kmers) - 1):
            if kmers[i + 1].startswith(kmers[i]):
                l.append(i+1)
        for i in reversed(range(len(l))):
            kmers.pop(l[i-1])
        """
        for i in kmers:
            e = edge(i)
            self.add_edge(e)

    def add_edge(self, edge):
        """ Create the edges and update the degrees for the nodes """
        if edge.label not in self.edges.keys():
            self.edges[edge.label] = edge
            self.nodes[edge.label[:-1]].edges.append(edge)
        else:
            self.edges[edge.label].multiplicity += 1

        self.nodes[edge.label[1:]].in_degree += 1
        self.nodes[edge.label[:-1]].out_degree += 1

    def labels_from_kmers(self, k, kmers):
        l = sorted(set([i[:k-1] for i in kmers] + [i[1:k] for i in kmers]))
        return l

    def export_graphviz(self):
        result = ''
        result += 'digraph {\n'
        result += '   graph [nodesep=2, size="300,300"];\n'
        for i in self.nodes:
            node = self.nodes[i]
            result += '    N%d [shape="box", style="rounded", label="%s"];\n' % (node.index, node.label)

        for i in self.edges:
            src = self.nodes[i[:-1]].index
            dst = self.nodes[i[1:]].index
            result += '    N%d -> N%d' % (src, dst)
            if (len(i) > 0):
                result += ' [label="%s", penwidth=1.0]' % i
            result += ';\n'
        result += '    overlap=false;\n'
        result += '}\n'
        return result


    def get_start(self):
        for n in self.nodes:
            if self.nodes[n].get_degree() == -1 and self.nodes[n].in_degree == 0:
                return self.nodes[n]

    def get_all_contigs(self, node = None):
        # Start at the beginning if we haven't started yet
        if node == None:
            node = self.get_start()

        node = self.get_contig(node)

        # If we've hit the end of the graph, start collapsing, else, continue
        count = 20
        while node != False and count > 0:
            count -= 1
            for edge in node.edges:
                node = self.get_contig(self.nodes[edge.label[:-1]],self.nodes[edge.label[1:]])

    def get_contig(self,start = None, next = None):
        node = start
        if next == None:
            path = []
        else:
            node = next
            path = [start]

        while node.out_degree <= 1:
            path.append(node)
            if (len(node.edges) > 0):
                node = self.nodes[node.edges[0].label[1:]]
            else:
                self.contigs[path[0].label] = path
                return False
        # Add the final node that broke the while condition
        path.append(node)
        self.contigs[path[0].label + path[1].label[-1]] = path
        return node

    def save_graph(self):
        file = open("output.gv","w")
        file.write(self.export_graphviz())
        file.close()
        print("[DONE]   ->  Saved Graph")

k = 21
kmers = kmer.create_kmers(k)

g = graph(k,kmers)
text = "fast. Nobody saw it "
text = "for a squid, fried e"
g.get_all_contigs()
g.save_graph()
print("[DONE]   ->  Contigs created")
