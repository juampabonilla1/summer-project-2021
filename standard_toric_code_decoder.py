import networkx as nx

import minimum_weight_perfect_matching as mwpm
from weight_functions import weight_of_edge_standard_lattice


def standard_decode(glued_lattice, p, bias):
    """
    Decodes biased Pauli noise on a standard square periodic toric code
    lattice.

    :param glued_lattice: A GluedLattices object which contains the standard
    toric code lattice as one of its attributes.
    :type glued_lattice: object
    :param p: The physical error probability on any one qubit.
    :type p: float
    :param bias: The bias coefficient of the noise toward dephasing.
    :type bias: float
    """
    L = glued_lattice.size

    # Match along even (white) and odd (black) parity plaquettes (in terms of
    # the sum of x and y indices) separately.
    path_weights = (1, 1, 1)
    nodes_w, nodes_b = [], []
    for i in range(L):
        for j in range(L):
            if glued_lattice.lattice_standard.plaquettes[i][j].state == 1:
                if (i + j) % 2 == 0:
                    nodes_w.append(str(i) + ',' + str(j))
                else:
                    nodes_b.append(str(i) + ',' + str(j))
    G = nx.Graph()
    for v1 in nodes_w:
        v1_1 = v1 + ',1'
        for v2 in nodes_w:
            v1_2 = v2 + ',1'
            if v1 != v2:
                wt_11_12 = weight_of_edge_standard_lattice(L, v1_1, v1_2,
                                                           path_weights)
                G.add_weighted_edges_from([(v1_1, v1_2, wt_11_12)])
    for v1 in nodes_b:
        v1_1 = v1 + ',1'
        for v2 in nodes_b:
            v1_2 = v2 + ',1'
            if v1 != v2:
                wt_11_12 = weight_of_edge_standard_lattice(L, v1_1, v1_2,
                                                           path_weights)
                G.add_weighted_edges_from([(v1_1, v1_2, wt_11_12)])

    matching = mwpm.mwpm(G)

    for u, v in matching:
        glued_lattice.lattice_standard.correct_in_lattice1(u, v, 'standard')
