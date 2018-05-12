class edge:
    def __init__(self,start,end):
        self.start = start
        self.end = end
        self.visited = False

# Given a string 'str' and a size of a k-mer 'k', this will construct the set of nodes/edges for a De Bruijn graph
class graph:
    def __init__(self, kmer):
        self.edges = None
        self.nodes = None
        self.kmer = kmer
        self.start_node = None
        self.end_node = None

    def create_graph(self, str):
        # If we're not appending to a graph, these will be empty
        k = self.kmer

        if self.edges == None:
            self.edges = []
        if self.nodes == None:
            self.nodes = set()

        for i in range(len(str) - k + 1):
            lkmer = str[i:i+k-1]
            rkmer = str[i+1:i+k]
            # Store the starting node
            if self.start_node == None:
                self.start_node = lkmer
                self.end_node = rkmer
            else:
                # update our start/end nodes
                if rkmer == self.start_node:
                    self.start_node = lkmer
                if lkmer == self.end_node:
                    self.end_node = rkmer
            # Creating left/right k-mers
            self.edges.append([str[i:i+k-1], str[i+1:i+k]])
            self.nodes.add(str[i:i+k-1])
            self.nodes.add(str[i+1:i+k])

        self.edges.append([self.end_node, self.start_node])


        return self.nodes, self.edges

    def update_graph(self, str):
        k = self.kmer
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

    def get_euler(self, graph):
        edges = graph
        graph = {}
        degree = {}
        start = edges[0][0]
        count_e = 0
        for e in edges:
            if not e[0] in graph:
                graph[e[0]] = {}
            if not e[0] in degree:
                degree[e[0]] = 0
            if not e[1] in graph:
                graph[e[1]] = {}
            if not e[1] in degree:
                degree[e[1]] = 0
            graph[e[0]][e[1]] = 1
            graph[e[1]][e[0]] = 1
            degree[e[0]] += 1
            degree[e[1]] += 1
            count_e += 1
        max_d = 0
        this_ = 0
        for v, d in degree.items():
            if not d%2 == 0 and d > 1:
                # Eulerian tour not possible as odd degree found!
                return False
            if d>max_d:
                this_ = v
                max_d = d
        visited_e = {}
        def is_visited(i, j):
            key = str(sorted([i,j]))
            if key in visited_e:
                return True
            else:
                visited_e[key] = True
                return False
        start = this_
        route = [start]
        indexof = {}
        indexof[start] = 0
        while count_e>0:
            flag = False
            for to_v in graph[this_]:
                if not is_visited(to_v, this_):
                    route.append([to_v])
                    indexof[to_v] = len(route)-1
                    degree[to_v] -= 1
                    if degree[to_v] == 0:
                        del degree[to_v]
                    degree[this_] -= 1
                    if degree[this_] == 0:
                        del degree[this_]
                    this_ = to_v
                    flag = True
                    count_e -= 1
                    break
            if not flag:
                break
        for key, v in degree.items():
            if v <=0:
                continue
            try:
                ind = indexof[key]
            except Exception as e:
                continue
            this_ = key
            while count_e>0:
                flag = False
                for to_v in graph[this_]:
                    if not is_visited(to_v, this_):
                        route[ind].append(to_v)
                        degree[to_v] -= 1
                        degree[this_] -= 1
                        this_ = to_v
                        flag = True
                        count_e -= 1
                        break
                if not flag:
                    break
        route_ref = []
        for r in route:
            if type(r) == list:
                for _r in r:
                    route_ref.append(_r)
            else:
                route_ref.append(r)
        return route_ref

    def concat_euler(self):
        print (self.get_euler())

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
