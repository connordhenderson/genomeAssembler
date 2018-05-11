import networkx as nx
import matplotlib.pyplot as plt

def draw(edges):
    G = nx.DiGraph()
    G.add_edges_from(edges)

    val_map = {'A': 1.0,
               'D': 0.5714285714285714,
               'H': 0.0}

    values = [val_map.get(node, 1.25) for node in G.nodes()]

    # Specify the edges you want here
    edge_colours = ['black' if not edge in red_edges else 'red'
                    for edge in G.edges()]
    black_edges = [edge for edge in G.edges() if edge not in red_edges]

    # Need to create a layout when doing
    # separate calls to draw nodes and edges
    pos = nx.spring_layout(G)
    plt.figure(1, figsize=(11,6.5))
    nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'),
                           node_color = values, node_size = 25)
    #nx.draw_networkx_labels(G, pos, font_size=6)
    nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=True)
    plt.show()
