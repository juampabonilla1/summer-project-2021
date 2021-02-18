import copy
import itertools

import SquareLattice1 as sl
from weight_functions import weight_of_edge
from graph_functions import perform_mwpm, get_weight_of_edge, create_graph, \
    add_weighted_edge_to_graph
from helper_functions import are_coords_equal, is_coord_in_tuple


def decode_glued_lattice(glued_lattice, path_weights, g_type='nx'):
    """
    Decode Pauli noise on a d * 2d lattice which is a concatenation of 2
    d * d lattices.

    :param glued_lattice: The concatenated lattice with Pauli noise on
    its qubits.
    :type glued_lattice: object
    :param path_weights: The weights of edges in the main and conditional
    panels. In the form path_weights=(A, B1, B2) where A, B1 and B2 are the
    weights of an edge in the main panel, low weight horizontal or vertical
    edge in the conditional panel and high weight edge in the conditional
    panel, respectively.
    :type path_weights: tuple of float
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools.
    :type g_type: str
    """
    L = glued_lattice.size

    # Build graph for matching.
    defects = create_list_of_defects_given_syndrome(glued_lattice)
    main_results = build_graph_from_list_of_defects(defects, L,
                                                    path_weights,
                                                    is_cut=False,
                                                    all_parity=False,
                                                    g_type=g_type)
    G, main_span, main_correc_defc = main_results
    # Match.
    matching = perform_mwpm(G, g_type)
    # matching_copy = copy.deepcopy(matching)

    # Create dummy lattice with cut.
    y_cut = int(L / 2)
    assert y_cut > 0 and y_cut < L  # Don't cut at the boundaries.

    # Remove defects whose path crosses cut.
    removed = remove_defects_whose_path_crosses_cut(matching, y_cut, L,
                                                    defects, path_weights)
    dummy_graph_defects, nodes_removed = removed
    nodes_removed_span = get_span_of_nodes_removed(nodes_removed, main_span)
    nds_rmvd_correc_defc = get_correc_defc_of_nds_rmvd(nodes_removed,
                                                       main_correc_defc)

    # Construct dummy graph.
    dummy_results = build_graph_from_list_of_defects(dummy_graph_defects, L,
                                                     path_weights,
                                                     is_cut=True,
                                                     y_cut=y_cut,
                                                     all_parity=True,
                                                     include_dummies=True,
                                                     g_type=g_type)
    dummy_G, dummy_span, dummy_correc_defc = dummy_results

    # Match on dummy graph.
    dummy_matching = perform_mwpm(dummy_G, g_type)

    # Finalise the inequivalent matching.
    inequiv_results = find_inequiv_matching(dummy_matching,
                                            dummy_G,
                                            nodes_removed,
                                            y_cut, L,
                                            path_weights,
                                            dummy_span,
                                            dummy_correc_defc,
                                            nodes_removed_span,
                                            nds_rmvd_correc_defc,
                                            g_type)
    inequiv_matching, inequiv_span, inequiv_correc_defc = inequiv_results
    # inequiv_matching_copy = copy.deepcopy(inequiv_matching)

    wt1, corr1 = create_sub_networks_and_find_weight(matching, main_span,
                                                     main_correc_defc, L,
                                                     g_type)
    wt2, corr2 = create_sub_networks_and_find_weight(inequiv_matching,
                                                     inequiv_span,
                                                     inequiv_correc_defc, L,
                                                     g_type)
    """
    print(matching_copy)
    print()
    print(inequiv_matching_copy)
    print()
    print(wt1, wt2)
    """

    if wt1 <= wt2:
        for correction in corr1:
            glued_lattice.lattice1.apply_correction_from_lst(correction)
    else:
        for correction in corr2:
            glued_lattice.lattice1.apply_correction_from_lst(correction)

    # for u, v in matching_copy:
    #    glued_lattice.lattice1.correct_in_lattice1(u, v)


def get_span_of_nodes_removed(nodes_removed, main_span):
    """
    Gets the x and y coordinates spanned by the path associated with the nodes
    removed. Nodes were removed from the main matching because their path
    crossed the cut.

    :param nodes_removed: List of pairs of nodes that were removed from main
    matching. In the form nodes_removed=[(u1, v1), (u2, v2), ...] where each
    u1='x1,y1,lat1', v1='x2,y2,lat2'. x1, y1, x2, y2 are the appropriate
    coordinates of the defects and lat1, lat2 are the lattice where they
    reside. 1 for main lattice, 2 for conditional lattice, 3 for dummy defect
    (in conditional lattice always).
    :type nodes_removed: lst of tuple of str
    :param main_span: Dictionary of the span associated with every path
    between matched defects in the main original matching. Keys are defect
    pairs, (u, v), values are the x and y span of the path between them, in
    the form [y_span, x_span]. e.g. y_span=[1,2,3].
    :type main_span: dict
    :return: Dictionary with keys as the pair of removed nodes and values
    equal to the span of x and y coordinates crossed by the path, in the form
    [y_span, x_span], e.g. y_span=[1,2,3].
    :rtype: dict
    """
    nodes_removed_span = {}
    for u, v in nodes_removed:
        span = main_span[(u, v)]
        nodes_removed_span.update({(u, v): span, (v, u): span})
    return nodes_removed_span


def get_correc_defc_of_nds_rmvd(nodes_removed, main_correc_defc):
    """
    Gets the associated defects to correct in the main lattice given two
    defects, which were removed because the path between them crossed the
    cut on the conditional lattice.

    Notes:

    * If defect 1 is in the main lattice and defect 2 is in the main lattice,
      the associated defects to correct are the same defect 1 and defect 2.
    * If defect 1 is in the main lattice and defect 2 is in the conditional
      lattice, the associated defects to correct are defect 1 and the
      appropriate coordinate of the intermediate boundary plaquette crossed
      by the path between the defects.
    * If defect 1 and defect 2 are in the conditional lattice, the associated
      defects are defect 1 and defect 2, which are ignored by the correction
      code anyway, since the correction code only looks at the restriction of
      the path on the main lattice.

    :param nodes_removed: List of pairs of nodes that were removed from main
    matching. In the form nodes_removed=[(u1, v1), (u2, v2), ...] where each
    u1='x1,y1,lat1', v1='x2,y2,lat2'. x1, y1, x2, y2 are the appropriate
    coordinates of the defects and lat1, lat2 are the lattice where they
    reside. 1 for main lattice, 2 for conditional lattice, 3 for dummy defect
    (in conditional lattice always).
    :type nodes_removed: lst of tuple of str
    :param main_correc_defc: Dictionary of the 'defects to correct' associated
    with every pair of defects in the main original matching. Keys are the
    matched defects, (u, v), values are the associated nodes to correct, see
    Notes above.
    :type main_correc_defc: dict
    :return: Dictionary with keys as the pair of removed nodes and values
    equal to the span of x and y coordinates crossed by the path, in the form
    [y_span, x_span], e.g. y_span=[1,2,3].
    :rtype: dict
    """
    nodes_removed_correc_defc = {}
    for u, v in nodes_removed:
        correc_defc = main_correc_defc[(u, v)]
        nodes_removed_correc_defc.update({(u, v): correc_defc,
                                          (v, u): correc_defc})
    return nodes_removed_correc_defc


def remove_defects_whose_path_crosses_cut(matching, y_cut, L,
                                          dummy_graph_defects,
                                          path_weights):
    """
    Removed the pair of defects from the main matching whose path crosses the
    cut in the conditional lattice.

    :param matching: The main original matching. In the form
    matching=[(u1, v1), (u2, v2), ...] where each u='x1,y1,lat1' and each
    v='x2,y2,lat2' are the defects. x1, y1, x2, y2 are the appropriate
    coordinates of the defects and lat1, lat2 are the lattice where they
    reside. 1 for main lattice, 2 for conditional lattice.
    :type matching: lst of tuple of str
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice.
    :type y_cut: int
    :param L: Size of LxL toric code lattice used in the concatenation.
    :type L: int
    :param dummy_graph_defects: List of defects which will be used for the
    inequivalent matching, assuming a cut on the conditional lattice. In the
    form [u1, u2, ...], each u='x,y,lat'.
    :type dummy_graph_defects: lst
    :param path_weights: The weights of edges in the main and conditional
    panels. In the form path_weights=(A, B1, B2) where A, B1 and B2 are the
    weights of an edge in the main panel, low weight horizontal or vertical
    edge in the conditional panel and high weight edge in the conditional
    panel, respectively.
    :type path_weights: tuple of float
    :return: The updated dummy graph defects after removing the defects
    associated with paths in the original matching which cross the cut, as
    well as a list of the removed defects, grouped in pairs as they appeared
    in the original matching.
    :rtype: tuple of lst and lst.
    """
    # List of ordered tuples. First coordinate in tuple is to left of cut.
    # Second coordinate is to right of cut.
    nodes_removed = []

    for u, v in matching:
        x1, y1, lat1 = v.split(",")
        x2, y2, lat2 = u.split(",")
        x1, y1, lat1 = int(x1), int(y1), int(lat1)
        x2, y2, lat2 = int(x2), int(y2), int(lat2)
        rmv = False
        if lat1 == lat2:
            if lat1 == 2:
                if (y1 < y_cut and y2 >= y_cut):
                    rmv = True
                    nodes_removed.append((v, u))
                elif (y1 >= y_cut and y2 < y_cut):
                    rmv = True
                    nodes_removed.append((u, v))
        else:
            # Make (x1, y1) the coordinates in lattice 1 (main).
            if lat1 == 2:
                x1, x2 = x2, x1
                y1, y2 = y2, y1
                u, v = v, u
            _, _, through_bd = weight_of_edge(L, v, u, path_weights,
                                              ret_interm_plaq=True)
            if (through_bd == 'middle' and y2 >= y_cut):
                rmv = True
                nodes_removed.append((v, u))
            elif (through_bd == 'edge' and y2 < y_cut):
                rmv = True
                nodes_removed.append((u, v))
        if rmv:
            dummy_graph_defects.remove(u)
            dummy_graph_defects.remove(v)

    return (dummy_graph_defects, nodes_removed)


def create_sub_networks_and_find_weight(matching, span, correc_defc, L,
                                        g_type):
    """
    Given a matching on the concatenated glued lattice, finds the networks
    of connected defects from the amalgamation of paths in the main and
    conditional lattice, and calculates the weight of each network as well
    as returns the minimal weight recovery operator to correct the network.

    :param matching: The main original matching. In the form
    matching=[(u1, v1), (u2, v2), ...] where each u='x1,y1,lat1' and each
    v='x2,y2,lat2' are the defects. x1, y1, x2, y2 are the appropriate
    coordinates of the defects and lat1, lat2 are the lattice where they
    reside. 1 for main lattice, 2 for conditional lattice.
    :type matching: lst of tuple of str
    :param span: Dictionary of the span associated with every path
    between matched defects. Keys are defect pairs, (u, v), values are the
    x and y span of the path between them, in the form [y_span, x_span].
    e.g. y_span=[1,2,3].
    :type span: dict
    :param correc_defc: Dictionary of the 'defects to correct' associated
    with every pair of defects in the matching. Keys are the matched defects,
    (u, v), values are the associated nodes to correct, see Notes in
    func:`get_correc_defc_of_nds_rmvd`.
    :type correc_defc: dict
    :param L: Size of LxL toric code lattice used in the concatenation.
    :type L: int
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools.
    :type g_type: str
    :return: The total weight of the matching as well as the recovery operator
    to correct every pair of defects in the matching. The total correction is
    in the form all_corrections=[corr1, corr2, ...] where each
    corr=[X_qubits, Y_qubits, Z_qubits] and each ?_qubits is a list of
    coordinates where to apply the corresponding Pauli operator.
    :rtype: tuple of int and lst
    """
    all_network_weights = []
    all_networks = []
    all_corrections = []

    while True:
        if len(matching) == 0:
            break

        network = []

        u, v = matching[0]

        network.append((u, v))
        matching.remove((u, v))

        if are_coords_equal(u, v):
            pass
        else:
            target, to_retrieve = u, v
            while True:
                pair = list(filter(lambda x: is_coord_in_tuple(to_retrieve, x),
                                   matching))
                assert len(pair) == 1
                u, v = pair[0]
                # Make v the defect to retrieve.
                if are_coords_equal(u, to_retrieve):
                    u, v = v, u
                else:
                    assert are_coords_equal(v, to_retrieve)

                network.append((u, v))
                matching.remove(pair[0])

                if are_coords_equal(u, target):
                    break

                to_retrieve = u

        all_networks.append(network)

        # print("network = ")
        # print(network)

        wt, correction = weight_and_correction_of_network(network, span,
                                                          correc_defc, L)

        """
        print("weight = ")
        print(wt)
        print("_______")
        """

        all_network_weights.append(wt)
        all_corrections.append(correction)
    return (sum(all_network_weights), all_corrections)


def weight_and_correction_of_network(network, span, correc_defc, L):
    """
    Calculates the weight of a single network and finds the recovery operator
    which corrects it.

    Notes:

    * This function should be changed depending on how the function to
      calculate the weights of the networks will be implemented.

    :param network: List of pairs of defects which make up the network. In the
    form network=[(u1, v1), (u2, v2), ...].
    :type network: lst
    :param span: Dictionary of the span associated with every path
    between matched defects. Keys are defect pairs, (u, v), values are the
    x and y span of the path between them, in the form [y_span, x_span].
    e.g. y_span=[1,2,3].
    :type span: dict
    :param correc_defc: Dictionary of the 'defects to correct' associated
    with every pair of defects in the matching. Keys are the matched defects,
    (u, v), values are the associated nodes to correct, see Notes in
    func:`get_correc_defc_of_nds_rmvd`.
    :type correc_defc: dict
    :param L: Size of LxL toric code lattice used in the concatenation.
    :type L: int
    :return: The weight of the network as well as the recovery operator
    to correct every pair of defects in the network, in the form
    final_correction=[X_qubits, Y_qubits, Z_qubits] and each ?_qubits is a
    list of coordinates where to apply the corresponding Pauli operator.
    :rtype: tuple of int and lst
    """
    # Find the edges of the network.
    all_y_sp, all_x_sp = set(), set()
    for u, v in network:
        y_sp, x_sp = span[(u, v)]
        all_y_sp |= set(y_sp)
        all_x_sp |= set(x_sp)
    all_y_sp, all_x_sp = list(all_y_sp), list(all_x_sp)

    if len(all_x_sp) > int(L / 2) or len(all_y_sp) > int(L / 2):
        # Create a lattice and apply the naive correction.
        lattice = sl.SquareLattice1(L)
        for u, v in network:
            u_to_correc, v_to_correc = correc_defc[(u, v)]
            # print("u, v = " + str((u, v)))
            # print("u_to_correct, v_to_correct = "
            #        + str((u_to_correc, v_to_correc)))
            lattice.correct_in_lattice1(u_to_correc, v_to_correc)
            # print(lattice)
            # lattice.print_plaquettes()

        # print(lattice)

        final_wt = lattice.weight()
        final_correction = lattice.cast_qubit_state_to_list()
        return (final_wt, final_correction)

    # Find the stabilisers within the network boundaries.
    all_encl_stabs = []
    for i in range(L):
        for j in range(L):
            left_bd = j in all_y_sp
            top_bd = i in all_x_sp
            right_bd = (j + 1) % L in all_y_sp
            bottom_bd = (i + 1) % L in all_x_sp
            if left_bd and top_bd and right_bd and bottom_bd:
                all_encl_stabs.append((i, j))

    # Find all 2^S combinations of on or off stabilisers from the S
    # stabilisers enclosed.
    all_stab_combs = []
    for tot_on_stabs in range(1, len(all_encl_stabs) + 1):
        all_stab_combs += list(itertools.combinations(all_encl_stabs,
                                                      tot_on_stabs))
    # Create a lattice and apply the naive correction.
    lattice = sl.SquareLattice1(L)
    for u, v in network:
        u_to_correc, v_to_correc = correc_defc[(u, v)]
        lattice.correct_in_lattice1(u_to_correc, v_to_correc)

    final_wt = lattice.weight()
    final_correction = lattice.cast_qubit_state_to_list()

    # Apply every possible combination of on or off stabilisers and find the
    # least weight correction.
    for on_stabs in all_stab_combs:
        trial_lattice = copy.deepcopy(lattice)
        for stab in on_stabs:
            trial_lattice.apply_stabiliser(stab)
        wt = trial_lattice.weight()
        if wt < final_wt:
            final_wt = wt
            final_correction = trial_lattice.cast_qubit_state_to_list()

    return (final_wt, final_correction)


def find_inequiv_matching(dummy_matching, dummy_G, nodes_removed, y_cut,
                          L, path_weights, dummy_span, dummy_correc_defc,
                          nodes_removed_span, nds_rmvd_correc_defc, g_type):
    """
    Finds the inequivalent matching around the lattice. Recovers the removed
    nodes and finalises the inequivalent/alternate matching.

    :param dummy_matching: Matching conducted on the defects of the lattice,
    assuming it was cut. In the form [(u1, v1), (u2, v2), ...]
    :type dummy_matching: lst
    :param dummy_G: The graph of defects and the weight between every pair,
    assuming a cut on the lattice. The distance between every defect and each
    dummy defect is the weighted horizontal distance to the cut.
    :type dummy_G: nx.Graph() or qecsim simpleGraph() or lst (see g_type).
    :param nodes_removed: List of pairs of nodes that were removed from main
    matching. In the form nodes_removed=[(u1, v1), (u2, v2), ...] where each
    u1='x1,y1,lat1', v1='x2,y2,lat2'.
    :type nodes_removed: lst of tuple of str
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice.
    :type y_cut: int
    :param L: Size of LxL toric code lattice used in the concatenation.
    :type L: int
    :param path_weights: The weights of edges in the main and conditional
    panels. In the form path_weights=(A, B1, B2) where A, B1 and B2 are the
    weights of an edge in the main panel, low weight horizontal or vertical
    edge in the conditional panel and high weight edge in the conditional
    panel, respectively.
    :type path_weights: tuple of float
    :param dummy_span: The x and y coordinates spanned by the paths between
    matched pairs. Dictionary with keys the matched pairs, (u, v), and values
    the y and x span, [y_span, x_span]. e.g. y_span=[1,2,3].
    :type dummy_span: dict
    :param dummy_correc_defc: The pair of defects to correct in lattice 1
    associated to every matched pair. See the note in
    the :func:`get_correc_defc_of_nds_rmvd`.
    :type dummy_correc_defc: dict
    :param nodes_removed_span: The y and x span of the nodes that were removed
    from the original matching because the path joining them crossed the cut.
    :type nodes_removed_span: dict
    :param nds_rmvd_correc_defc: The defects to correct associated to the
    removed nodes.
    :type nds_rmvd_correc_defc: dict
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools, 'blsm' for a
    direct used of Blossom V.
    :type g_type: str
    :return: The inequivalent matching as a list of pairs of defects,
    [(u1, v1), (u2, v2), ...], the y and x span of the paths between matched
    defects as a dict with keys the pair, (u, v), and values [y_span, x_span],
    and the pair of defects to correct in lattice 1 associated to every
    matched pair. as a dict (See the note in the
    :func:`get_correc_defc_of_nds_rmvd`).
    :rtype: tuple of lst, dict, dict
    """
    inequivalent_matching = []
    inequivalent_span = {}
    inequiv_correc_defc = {}

    # Find the defects the dummies matched to, remove the dummies and join
    # these remnant defects
    dm1_matching = list(filter(lambda x: x[0] == 'dm1' or x[1] == 'dm1',
                        dummy_matching))[0]
    if dm1_matching not in [('dm1', 'dm2'), ('dm2', 'dm1')]:
        v_to_dm1 = get_partner_to_dm(dm1_matching)
        wt1 = get_weight_of_edge(dummy_G, dm1_matching[0], dm1_matching[1],
                                 g_type)

        dm2_matching = list(filter(lambda x: x[0] == 'dm2' or x[1] == 'dm2',
                            dummy_matching))[0]
        v_to_dm2 = get_partner_to_dm(dm2_matching)
        wt2 = get_weight_of_edge(dummy_G, dm2_matching[0], dm2_matching[1],
                                 g_type)

        dummy_matching.remove(dm1_matching)
        dummy_matching.remove(dm2_matching)

        wt_btwn_dummies = weight_between_dummies_given_their_match(v_to_dm1,
                                                                   v_to_dm2,
                                                                   L, y_cut,
                                                                   path_weights
                                                                   )
        wt_btwn_remnant_defects = wt1 + wt2 + wt_btwn_dummies

        # Go through the list of removed defects and see if joining these
        # defects to the remnant defects is lower weight than joining remnant
        # defects together.
        left_match, right_match = None, None
        for u, v in nodes_removed:
            wt_to_v1 = weight_of_edge(L, u, v_to_dm1, path_weights,
                                      is_cut=True, y_cut=y_cut)
            wt_to_v2 = weight_of_edge(L, v, v_to_dm2, path_weights,
                                      is_cut=True, y_cut=y_cut)
            total_wt = wt_to_v1 + wt_to_v2
            # NOTE: potential for problems below. Should it be:
            # + weight of the removed path, that is, between u and v?
            if total_wt <= wt_btwn_remnant_defects:
                left_match, right_match = u, v
                wt_btwn_remnant_defects = total_wt
        if left_match is not None:
            assert right_match is not None

            nodes_removed.remove((left_match, right_match))
            nodes_removed_span.pop((left_match, right_match), None)
            nodes_removed_span.pop((right_match, left_match), None)
            nds_rmvd_correc_defc.pop((left_match, right_match), None)
            nds_rmvd_correc_defc.pop((right_match, left_match), None)

            inequivalent_matching.append((left_match, v_to_dm1))
            inequivalent_matching.append((right_match, v_to_dm2))

            update_span_with_final_match(inequivalent_span, left_match,
                                         v_to_dm1, L, y_cut, path_weights,
                                         is_cut=True)
            update_span_with_final_match(inequivalent_span, right_match,
                                         v_to_dm2, L, y_cut, path_weights,
                                         is_cut=True)

            update_correc_defc_with_final_match(inequiv_correc_defc,
                                                left_match, v_to_dm1, L,
                                                y_cut, path_weights,
                                                is_cut=True)
            update_correc_defc_with_final_match(inequiv_correc_defc,
                                                right_match, v_to_dm2, L,
                                                y_cut, path_weights,
                                                is_cut=True)

        else:
            dm1, dm2 = coord_of_dummies_given_their_match(v_to_dm1, v_to_dm2,
                                                          L, y_cut)
            inequivalent_matching.append((v_to_dm1, dm1))
            inequivalent_matching.append((dm1, dm2))
            inequivalent_matching.append((v_to_dm2, dm2))

            update_span_with_final_match(inequivalent_span, v_to_dm1,
                                         dm1, L, y_cut, path_weights,
                                         is_cut=True)
            update_span_with_final_match(inequivalent_span, dm1,
                                         dm2, L, y_cut, path_weights,
                                         is_cut=False)
            update_span_with_final_match(inequivalent_span, v_to_dm2,
                                         dm2, L, y_cut, path_weights,
                                         is_cut=True)

            update_correc_defc_with_final_match(inequiv_correc_defc,
                                                v_to_dm1, dm1, L,
                                                y_cut, path_weights,
                                                is_cut=True)
            update_correc_defc_with_final_match(inequiv_correc_defc,
                                                dm1, dm2, L,
                                                y_cut, path_weights,
                                                is_cut=True)
            update_correc_defc_with_final_match(inequiv_correc_defc,
                                                v_to_dm2, dm2, L,
                                                y_cut, path_weights,
                                                is_cut=True)
    else:
        dummy_matching.remove(dm1_matching)

    inequivalent_matching += dummy_matching + nodes_removed
    inequivalent_span.update(dummy_span)
    inequivalent_span.update(nodes_removed_span)
    inequiv_correc_defc.update(dummy_correc_defc)
    inequiv_correc_defc.update(nds_rmvd_correc_defc)
    return (inequivalent_matching, inequivalent_span, inequiv_correc_defc)


def update_span_with_final_match(inequivalent_span, u, v, L, y_cut,
                                 path_weights, is_cut):
    """
    Updates the dictionary of spanned defects by the paths joining the
    matched pairs.

    :param inequivalent_span: The y and x span of the paths between matched
    defects as a dict with keys the pair, (u, v), and values [y_span, x_span].
    :type inequivalent_span: dict
    :param u: The first defect, in the form 'x1,y1,lat1'.
    :type u: str
    :param v: The second defect, in the form 'x2,y2,lat2'.
    :type v: str
    :param L: The size of the LxL toric code lattice concatenated.
    :type L: int
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice, if the lattice was cut.
    :type y_cut: int
    :param path_weights: The weights of edges in the main and conditional
    panels.
    :type path_weights: tuple of float
    :param is_cut: Whether the conditional lattice is cut or not.
    :type is_cut: bool
    """
    _, y_sp, x_sp = weight_of_edge(L, u, v, path_weights, is_cut=is_cut,
                                   y_cut=y_cut, ret_span=True)
    inequivalent_span.update({(u, v): [y_sp, x_sp], (v, u): [y_sp, x_sp]})


def update_correc_defc_with_final_match(inequiv_correc_defc, u, v, L, y_cut,
                                        path_weights, is_cut):
    """
    Updates the dictionary of defects to correct in lattice 1
    associated to every matched pair in the inequivalent matching.

    :param inequivalent_span: The pair of defects to correct in lattice 1
    associated to every matched pair. See the note in
    the :func:`get_correc_defc_of_nds_rmvd`.
    :type inequivalent_span: dict
    :param u: The first defect, in the form 'x1,y1,lat1'.
    :type u: str
    :param v: The second defect, in the form 'x2,y2,lat2'.
    :type v: str
    :param L: The size of the LxL toric code lattice concatenated.
    :type L: int
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice, if the lattice was cut.
    :type y_cut: int
    :param path_weights: The weights of edges in the main and conditional
    panels.
    :type path_weights: tuple of float
    :param is_cut: Whether the conditional lattice is cut or not.
    :type is_cut: bool
    """
    _ = weight_of_edge(L, u, v, path_weights, is_cut=is_cut,
                       y_cut=y_cut, ret_span=True,
                       cor_defc=inequiv_correc_defc)


def weight_between_dummies_given_their_match(u, v, L, y_cut, path_weights):
    """
    Calculates the weighted distance between the dummy defects given the
    defects they matched to.

    :param u: Defect first dummy matched to, in the form 'x1,y1,lat1'.
    :type u: str
    :param v: Defect second dummy matched to, in the form 'x1,y1,lat1'.
    :type v: str
    :param L: The size of the LxL toric code lattice concatenated.
    :type L: int
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice, if the lattice was cut.
    :type y_cut: int
    :param path_weights: The weights of edges in the main and conditional
    panels.
    :type path_weights: tuple of float
    :return: The weighted distance between the dummy defects.
    :rtype: float
    """
    x1, x2 = int(u.split(",")[0]), int(v.split(",")[0])
    dm1 = str(x1) + ',' + str(y_cut - 1) + ',2'
    dm2 = str(x2) + ',' + str(y_cut) + ',2'
    return weight_of_edge(L, dm1, dm2, path_weights)


def coord_of_dummies_given_their_match(u, v, L, y_cut):
    """
    Finds the coordinates of the dummy defects in the second lattice, given
    the defects they matched to.

    :param u: Defect first dummy matched to, in the form 'x1,y1,lat1'.
    :type u: str
    :param v: Defect second dummy matched to, in the form 'x1,y1,lat1'.
    :type v: str
    :param L: The size of the LxL toric code lattice concatenated.
    :type L: int
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice, if the lattice was cut.
    :type y_cut: int
    :return: The coordinates of the dummies, in the form (dm1, dm2), where
    each dm='x,y,3'. The 3 flags that it is a dummy defect.
    :rtype: tuple of str, str
    """
    x1, x2 = int(u.split(",")[0]), int(v.split(",")[0])
    dm1 = str(x1) + ',' + str(y_cut - 1) + ',3'
    dm2 = str(x2) + ',' + str(y_cut) + ',3'
    return (dm1, dm2)


def get_partner_to_dm(dm_matching):
    """
    Retrieves the defect the dummy matched to.

    :param dm_matching: The tuple containing the dummy and the defect the
    the dummy matched to.
    :type dm_matching: tuple
    :return: The defect the dummy defect matched to.
    :rtype: str
    """
    u, v = dm_matching
    if u not in ['dm1', 'dm2']:
        assert v in ['dm1', 'dm2']
        return u
    else:
        return v


def create_list_of_defects_given_syndrome(lattice):
    """
    Collects the coordinates of the plaquette locations where there is a
    defect in the lattice.

    :param lattice: A SquareLattice1 object.
    :type lattice: object
    :return: A list of plaquette locations in the concatenated lattice object
    where there are defects. In the form [df1, df2, ...], where each
    df='x,y,lat'. x is the row number, y is the column number, lat is the
    lattice: 1 for the main lattice, 2 for the conditional lattice.
    :rtype: lst
    """
    L = lattice.size
    defects = []
    for i in range(L):
        for j in range(L):
            if lattice.lattice1.plaquettes[i][j].state == 1:
                defects.append(str(i) + ',' + str(j) + ',1')
            if lattice.lattice2.plaquettes[i][j].state == 1:
                defects.append(str(i) + ',' + str(j) + ',2')
    return defects


def build_graph_from_list_of_defects(defects, L, path_weights, is_cut=False,
                                     y_cut=None, all_parity=False,
                                     include_dummies=False, g_type='nx'):
    """
    Creates a graph object from a list of defect coordinates in the lattice.
    Nodes are the defect locations and weighted edges between nodes are the
    weighted distance between defects in the lattice.

    :param defects: List of defects in the form [u1, u2, ...] where each
    u='x,y,lat'.
    :type defects: lst
    :param L: The size of the LxL toric code lattice concatenated.
    :type L: int
    :param path_weights: The weights of edges in the main and conditional
    panels. In the form path_weights=(A, B1, B2) where A, B1 and B2 are the
    weights of an edge in the main panel, low weight horizontal or vertical
    edge in the conditional panel and high weight edge in the conditional
    panel, respectively.
    :type path_weights: tuple of float
    :param is_cut: Whether the lattice is cut or not.
    :type is_cut: bool
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice, if the lattice was cut.
    :type y_cut: int
    :param all_parity: Whether we allow matching between different parity
    defects (dark and light) when one defects lies in lattice 1 and the other
    in lattice 2. Note that we always (never) allow different parity matching
    when the two defects are both in lattice 2 (1).
    :param include_dummies: Whether we are including the dummy defects on
    either side of the cut or not.
    :type include_dummies: bool
    :param g_type: The type of graph and matching tools being used. 'nx' for
    networkx, 'qecsim' for qecsim graph and matching tools, 'blsm' for a
    direct used of Blossom V.
    :type g_type: str
    :return: The graph object containing the defects as nodes, a dictionary
    of the spanned coordinates by the paths joining each pair of defects, with
    keys the defects, (u, v), and values [y_span, x_span] where e.g.
    y_span=[1,2,3], as well as a dictionary of he pair of defects to correct
    in lattice 1 associated to every matched pair. See the note in
    the :func:`get_correc_defc_of_nds_rmvd`.
    :rtype: tuple
    """
    A, B1, _ = path_weights

    spanned_coords = {}
    cor_defc = {}
    G = create_graph(g_type)

    for v1 in defects:
        x1, y1, lat1 = v1.split(",")
        x1, y1, lat1 = int(x1), int(y1), int(lat1)
        for v2 in defects:
            x2, y2, lat2 = v2.split(",")
            x2, y2, lat2 = int(x2), int(y2), int(lat2)
            if v1 != v2:
                if lat1 == lat2:
                    if lat1 == 1:
                        # Can only match along same parity plaquettes.
                        if (x1 + y1) % 2 == (x2 + y2) % 2:
                            wt, y_sp, x_sp = weight_of_edge(L, v1, v2,
                                                            path_weights,
                                                            is_cut=is_cut,
                                                            y_cut=y_cut,
                                                            ret_span=True,
                                                            cor_defc=cor_defc)
                            add_weighted_edge_to_graph(G, v1, v2, wt, g_type)
                            spanned_coords.update({(v1, v2): [y_sp, x_sp],
                                                   (v2, v1): [y_sp, x_sp]})
                    else:
                        wt, y_sp, x_sp = weight_of_edge(L, v1, v2,
                                                        path_weights,
                                                        is_cut=is_cut,
                                                        y_cut=y_cut,
                                                        ret_span=True,
                                                        cor_defc=cor_defc)
                        add_weighted_edge_to_graph(G, v1, v2, wt, g_type)
                        spanned_coords.update({(v1, v2): [y_sp, x_sp],
                                               (v2, v1): [y_sp, x_sp]})
                else:
                    # NOTE: this assumes that we are allowing different parity
                    # matching across lattices.
                    if ((x1 + y1) % 2 == (x2 + y2) % 2) or all_parity:
                        wt, y_sp, x_sp = weight_of_edge(L, v1, v2,
                                                        path_weights,
                                                        is_cut=is_cut,
                                                        y_cut=y_cut,
                                                        ret_span=True,
                                                        cor_defc=cor_defc)
                        add_weighted_edge_to_graph(G, v1, v2, wt, g_type)
                        spanned_coords.update({(v1, v2): [y_sp, x_sp],
                                               (v2, v1): [y_sp, x_sp]})

    if include_dummies:
        add_weighted_edge_to_graph(G, 'dm1', 'dm2', A * L + B1 * L, g_type)
        for v in defects:
            wt1, wt2 = weights_to_dummies(v, L, y_cut, path_weights)
            add_weighted_edge_to_graph(G, 'dm1', v, wt1, g_type)
            add_weighted_edge_to_graph(G, 'dm2', v, wt2, g_type)
    return (G, spanned_coords, cor_defc)


def weights_to_dummies(v, L, y_cut, path_weights):
    """
    Finds the weight of a defect to each of the dummy defects.

    Notes:

    * This weight is currently given as the weighted horizontal distance from
      the defect to either side of the cut on the conditional lattice.

    :param v: The defect under consideration, in the form 'x,y,lat'.
    :type v: str
    :param L: The size of the LxL toric code lattice concatenated.
    :type L: int
    :param y_cut: The vertical coordinate of the cut in the conditional
    lattice, if the lattice was cut.
    :type y_cut: int
    :param path_weights: The weights of edges in the main and conditional
    panels.
    :type path_weights: tuple of float
    :return: The weighted distance to the dummy on the left of the cut, wt1,
    and the weighted distance to the dummy on the right of the cut, wt2.
    :rtype: tuple
    """
    A, B1, _ = path_weights

    x, y, lat = v.split(",")
    x, y, lat = int(x), int(y), int(lat)

    if lat == 1:
        wt1 = A * (L - y) + B1 * (y_cut - 1)
        wt2 = A * (y + 1) + B1 * (L - y_cut - 1)
    else:
        if y < y_cut:
            wt1 = B1 * (y_cut - y - 1)
            wt2 = L * A + B1 * L - wt1 - 1
        else:
            wt2 = B1 * (y - y_cut)
            wt1 = L * A + B1 * L - wt2 - 1
    return (wt1, wt2)
