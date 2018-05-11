# Given a string 'str' and a size of a k-mer 'k', this will construct the set of nodes/edges for a De Bruijn graph
def create_graph(nodes, edges, str,k):
    # If we're not appending to a graph, these will be empty
    if edges == None:
        edges = []
    if nodes == None:
        nodes = set()

    for i in range(len(str) - k + 1):
        # Creating left/right k-mers
        edges.append([str[i:i+k-1], str[i+1:i+k]])
        nodes.add(str[i:i+k-1])
        nodes.add(str[i+1:i+k])
    return nodes, edges

def update_graph(nodes, edges, str, k):
        for i in range(len(str) - k + 1):
            # Creating left/right k-mers
            edges.append([str[i:i+k-1], str[i+1:i+k]])
            nodes.add(str[i:i+k-1])
            nodes.add(str[i+1:i+k])
        return nodes, edges


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

# We're going to attempt to implement Hierholzer's algorithm
def get_eulerian(edges):
    end = False

    d = create_dict(edges)

    # we need to create a dictionary of our edges
    for e in edges:
        d[e[0]] = e[1]

    v = list(d)[0]
    stack = [v]
    circuit = []

    current = 0
    next = 0

    # count the amount of edges on each vertex
    e_count = dict()
    for vx in range(len(list(d))):
        e_count[vx] = len(d[vx])

    while len(c_path) >= 0:
        # if we have an exit edge, carry on
        if len(d[vx]) > 0:
            stack.append(current)

            next = d[current][0]
            # remove that edge
            arr = d[current]
            arr = arr[1:]
            d[current] = arr

            current = next
        else:
            circuit.append(current)
            current = stack[-1]
            stack.pop()
    return circuit
