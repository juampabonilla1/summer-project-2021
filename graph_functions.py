import networkx as nx
import matplotlib.pyplot as plt

import minimum_weight_perfect_matching as mwpm
# import pypm_example as pe
# import qecsim.graphtools as gt


def create_graph(g_type):
    """
    Creates a graph object.

    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools, 'blsm' for a
    direct used of Blossom V.
    :type g_type: str
    :return: The graph object initialised with no nodes.
    :rtype: nx.Graph() or qecsim simpleGraph() or lst
    """
    if g_type == 'nx':
        return nx.Graph()
    elif g_type == 'blsm':
        return []
    elif g_type == 'qecsim':
        return gt.SimpleGraph()


def get_weight_of_edge(G, u, v, g_type):
    """
    Retrieves the weight of an edge between 2 particular nodes in a graph
    object.

    :param G: The graph object. Its composition depends on which structure we
    use. See :param:`g_type`.
    :type G: nx.Graph() or qecsim simpleGraph() or lst.
    :param u: First node. In the form 'x1,y1,lat1'.
    :type u: str
    :param v: Second node. In the form 'x2,y2,lat2'.
    :type v: str
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools, 'blsm' for a
    direct used of Blossom V.
    :type g_type: str
    """
    if g_type == 'nx':
        return G.edges[(u, v)]['weight']
    elif g_type == 'blsm':
        return list(filter(lambda x: (x[0] == u and x[1] == v) or
                                     (x[0] == v and x[1] == u), G))[0][2]
    elif g_type == 'qecsim':
        if (u, v) in G:
            return G[(u, v)]
        else:
            return G[(v, u)]


def perform_mwpm(G, g_type):
    """
    Perform minimum weight perfect matching on a weighted graph object.

    :param G: The graph object. Its composition depends on which structure we
    use. See :param:`g_type`.
    :type G: nx.Graph() or qecsim simpleGraph() or lst.
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools, 'blsm' for a
    direct used of Blossom V.
    :type g_type: str
    :return: The minimum weight perfect match in the form
    match=[(u1, v1), (u2, v2), ...] where each u, v are nodes in the graph
    corresponding to defects in the lattice.
    :rtype: lst
    """
    if g_type == 'nx':
        return mwpm.mwpm(G)
    elif g_type == 'blsm':
        return list(pe.mwpm(G))
    elif g_type == 'qecsim':
        return list(gt.mwpm(G))


def print_nx_graph(G):
    """
    Prints a networkx graph, including the weights of the connecting edges.

    :param G: The graph under consideration.
    :type G: nx.Graph()
    """
    plt.figure()
    pos = nx.spring_layout(G)
    nx.draw_networkx(G, pos)
    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
    plt.show()


def add_weighted_edge_to_graph(G, u, v, wt, g_type):
    """
    Adds two nodes to a graph and an edge between them with a predetermined
    weight.

    :param G: Graph object.
    :type G: nx.Graph() or qecsim simpleGraph() or lst
    :param u: First node.
    :type u: str
    :param v: Second node.
    :type v: str
    :param wt: Weight of edge.
    :type wt: float
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools, 'blsm' for a
    direct used of Blossom V.
    :type g_type: str
    """
    if g_type == 'nx':
        G.add_weighted_edges_from([(u, v, wt)])
    elif g_type == 'blsm':
        G.append((u, v, wt))
    elif g_type == 'qecsim':
        G.add_edge(u, v, wt)
