import kmer

class vertex:
    def __init__(self, label, index):
        self.index = index
        self.label = label
        self.outs = []

class edge:
    def __init__(self, l, r):
        self.left = l
        self.right = r
        self.label = ""

class graph:
    def __init__(self, kmers, k=9):
        """
        'k' should always be odd so that we don't return a compliment that
        is identical
        """
        if (k%2==0):
            k +=1
        vlist = nodes_from_kmers(kmers,k)

        self.k = k
        self.edges = []
        self.edge_labels = []
        self.indices = {v:i for i,v in enumerate(vlist)}
        self.vertices = {i:v for i,v in enumerate(vlist)}

        self.start = None
        self.end = None

        self.create_edges(kmers)

    """ Add a vertex by name (label) to the graph"""
    def add_vertex(self, label):
        i = len(self.indices)
        self.indices[label] = i
        self.vertices[i] = label

    """ Creates a directed edge between specified vertices"""
    def add_edge(self, vsrc, vdst, label='', repeats=False):
        e = (self.indices[vsrc], self.indices[vdst])
        if (repeats) or (e not in self.edges):
            self.edges.append(e)
            self.edge_labels.append(label)

    def create_edges(self, kmers):
        for i in kmers:
            self.add_edge(i[:self.k-1], i[1:self.k], i)

    """
    Gets the in/out degree of the nodes in the graph
    """
    def get_degrees(self):
        in_degree = {}
        out_degree = {}

        for src,dst in self.edges:
            out_degree[src] = out_degree.get(src, 0) + 1
            in_degree[dst] = in_degree.get(dst, 0) + 1

        return in_degree, out_degree

    def get_edges(self):
        return self.edges

    """
    Returns the Eulerian path (if one was found, False otherwise) as a list of
    edge ID's. Implemented using University of North Carolina's computational
    systems biology method
    """
    def get_eulerian_path(self):
        g = [(src,dst) for src,dst in self.edges]
        c_vertex = self.get_start()
        path = [c_vertex]
        # "next" is where vertices get inserted into our tour
        # it starts at the end (i.e. it is the same as appending),
        # but later "side-trips" will insert in the middle
        next = 1
        while len(g) > 0:
            for edge in g:
                if (edge[0] == c_vertex):

                    c_vertex = edge[1]
                    g.remove(edge)
                    path.insert(next, c_vertex)
                    next += 1
                    break
            else:
                for edge in g:
                    try:
                        next = path.index(edge[0]) + 1
                        c_vertex = edge[0]
                        break
                    except ValueError:
                        continue
                else:
                    print ("There is no path!")
                    return False
        return path

    """
    Given a list of edge IDs, this will return the corresponding list of edges.
    Implemented using the University of North Carolina's computational systems
    biology method
    """
    def get_labels_from_edge(self, path):
        edge_id = {}
        for i in range(len(self.edges)):
            edge_id[self.edges[i]] = edge_id.get(self.edges[i], []) + [i]
        edge_list = []
        for i in range(len(path) - 1):
            try:
                edge_list.append(self.edge_labels[edge_id[path[i],path[i+1]].pop()])
            except KeyError:
                return edge_list
        return edge_list

    def get_nodes(self):
        return self.vertices

    """
    Returns the starting node in the graph by comparing the amount of in/out
    degrees within the graph. This method will likely become unreliable as we
    start to get read errors, and should be reevaluated as a source of error
    """
    def get_start(self):
        in_degree, out_degree = self.get_degrees()
        start = -1
        end = -1

        for vert in self.vertices.keys():
            ins = in_degree.get(vert, 0)
            outs = out_degree.get(vert, 0)
            if (ins == outs):
                continue
            elif (ins - outs) == 1:
                end = vert
            elif (outs - ins) == 1:
                start = vert

        if (start >= 0) and (end >= 0):
            self.end = end
            self.start = start
            return start
        else:
            return -1

    def get_vertices(self):
        return self.vertices

    """
    Returns the graph as a string that is readable by GraphViz for visualizing and
    identifying sources of error for smaller graphs. Larger graphs can be created,
    but attempting to graph them will result in the graphing application to crash
    """
    def render(self, highlightPath = []):
        edge_id = {}
        for i in range(len(self.edges)):
            edge_id[self.edges[i]] = edge_id.get(self.edges[i], []) + [i]
        edge_set = set()

        for i in range(len(highlightPath) - 1):
            src = self.indices[highlightPath[i]]
            dst = self.indices[highlightPath[i+1]]
            edge_set.add(edge_id[src,dst].pop())
        result = ''
        result += 'digraph {\n'
        result += '   graph [nodesep=2, size="300,300"];\n'

        for index, label in self.vertices.items():
            result += '    N%d [shape="box", style="rounded", label="%s"];\n' % (index, label)
        for i, e in enumerate(self.edges):
            src, dst = e
            result += '    N%d -> N%d' % (src, dst)
            label = self.edge_labels[i]
            if (len(label) > 0):
                if (i in edge_set):
                    result += ' [label="%s", penwidth=3.0]' % (label)
                else:
                    result += ' [label="%s"]' % (label)
            elif (i in edge_set):
                result += ' [penwidth=3.0]'
            result += ';\n'
        result += '    overlap=false;\n'
        result += '}\n'
        return result

    """
    Gets the GraphViz readable string and saves it to a file labeled
    'output.txt'
    """
    def save_graph(self):
        file = open("output.gv","w")
        file.write(self.render())
        file.close()
        print("[DONE] Saved Graph")

    """
    A shortcut helper function that will get the Eulerian path and display it.
    Used when you don't need to explicitly print the edges/path (like debugging)
    """
    def print_euler(self):
        path = self.get_eulerian_path()
        path = self.get_labels_from_edge(path)
        seq = path[0][0:self.k]
        for kmer in path[1:]:
            seq += kmer[self.k-1]
        print (seq)

    """
    Gets the Eulerian tour and saves its sequence to a file, returns the sequence
    """
    def save_euler(self, fpath="Data/output.dat"):
        """Get the Eulerian tour, and concat the last letter of each node
        to give us the resulting sequence"""
        path = self.get_eulerian_path()
        if (path == False):
            return False
        else:
            print("Eulerian Path Found")

        path = self.get_labels_from_edge(path)
        seq = path[0][0:self.k]

        for kmer in path[1:]:
            seq += kmer[self.k-1]

        with open(fpath, 'w') as file:
            file.write(seq)

        return seq

"""
Used to construct the graph; ultimately redundant when you can just use the
constructor, but is a hold-over from previous methods that sees testing use
"""
def create_graph(kmers, k):
    g = None
    g = graph(kmers, k)
    print("[DONE]   ->  graph created")
    return g


"""
Creates a list of nodes from a list of k-mers
"""
def nodes_from_kmers(kmers, k):
    return sorted(set([i[:k-1] for i in kmers] + [i[1:k] for i in kmers]))
