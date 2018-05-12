class edge:
    def __init__(self,start,end):
        self.start = start
        self.end = end
        self.visited = False

# Given a string 'str' and a size of a k-mer 'k', this will construct the set of nodes/edges for a De Bruijn graph
class graph:
    def __init__(self):
        self.edges = None
        self.nodes = None

    def create_graph(self, str,k):
        # If we're not appending to a graph, these will be empty
        if self.edges == None:
            self.edges = []
        if self.nodes == None:
            self.nodes = set()

        for i in range(len(str) - k + 1):
            # Creating left/right k-mers
            self.edges.append([str[i:i+k-1], str[i+1:i+k]])
            self.nodes.add(str[i:i+k-1])
            self.nodes.add(str[i+1:i+k])
        return self.nodes, self.edges

    def update_graph(self, str, k):
            for i in range(len(str) - k + 1):
                # Creating left/right k-mers
                self.edges.append([str[i:i+k-1], str[i+1:i+k]])
                self.nodes.add(str[i:i+k-1])
                self.nodes.add(str[i+1:i+k])
            return self.nodes, self.edges

    # We're going to attempt to implement Hierholzer's algorithm
    def next_node(self, edge, current):
        return edge[0] if current == edge[1] else edge[1]

    def remove_edge(self, raw_list, discard):
        return [item for item in raw_list if item != discard]

    def get_eulerian(self):
        graph = self.edges
        search = [[[], graph[0][0], graph]]
        while search:
            path, node, unexplore = search.pop()
            path += [node]

            if not unexplore:
                return path

            for edge in unexplore:
                if node in edge:
                    search += [[path, self.next_node(edge, node), self.remove_edge(unexplore,edge)]]




# We need to be able to create a dictionary that will store information
def create_dict(edges):
    d = dict()
    for e in edges:
        if (e[0] in d):
            arr = []
            arr.append(d[e[0]])
            arr.append(e[1])
            d[e[0]] = arr
        else:
            d[e[0]] = [e[1]]
    return d
