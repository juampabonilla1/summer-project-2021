from operator import itemgetter


def weight_of_edge_standard_lattice(L, v, u, path_weights):
    """
    Calculates the weight of an edge between two defects in the standard
    square periodic toric code lattice.

    Notes:

    * This function assumes that we are decoding with a single square,
    periodic lattice. In particular, the path between two defects can go
    around the lattice. There is no concatenated conditional lattice next
    to it.

    :param L: Size of LxL periodic square lattice.
    :type L: int
    :param v: First defect. In the form 'x1,y1,lat1'.
    :type v: str
    :param u: Second defect. In the form 'x2,y2,lat2'.
    :type u: str
    :param path_weights: The weights of edges in the main and conditional
    panels. In the form path_weights=(A, B1, B2) where A, B1 and B2 are the
    weights of an edge in the main panel, low weight horizontal or vertical
    edge in the conditional panel and high weight edge in the conditional
    panel, respectively. In this function we only use A.
    :type path_weights: tuple of float
    :return: The weight of the edge in the lattice.
    :rtype: float
    """
    A, _, _ = path_weights

    x1, y1, lat1 = v.split(",")
    x2, y2, lat2 = u.split(",")
    x1, y1, lat1 = int(x1), int(y1), int(lat1)
    x2, y2, lat2 = int(x2), int(y2), int(lat2)

    assert lat1 == lat2 and lat1 == 1
    assert (x1 + y1) % 2 == (x2 + y2) % 2
    lx = min([abs(x2 - x1), L - abs(x2 - x1)])
    ly = min([abs(y2 - y1), L - abs(y2 - y1)])
    return A * max([lx, ly])


def weight_of_edge_main_lattice(L, x1, y1, x2, y2, A, ret_span=False):
    """
    Calculates the weight of an edge between two defects in the main lattice
    of a concatenated glued lattice object.

    :param L: Size of LxL toric code lattice used in the concatenation.
    :type L: int
    :param x1: Row coordinate of first defect.
    :type x1: int
    :param y1: Column coordinate of first defect.
    :type y1: int
    :param x2: Row coordinate of second defect.
    :type x2: int
    :param y2: Column coordinate of second defect.
    :type y2: int
    :param A: Weight of an edge in the main lattice.
    :type A: float
    :param ret_span: Whether to return the y and x coordinates spanned by the
    path between defects or not.
    :type ret_span: bool
    :return: The weight between defects and possibly the span of the path
    between them.
    :rtype: float or tuple
    """
    # In lattice1, weighted distance between 2 defects is the taxi
    # cab distance times weight of each step, A.

    # Make y1 the left coordinate.
    if y1 > y2:
        y1, y2 = y2, y1

    ly = abs(y2 - y1)
    y_span = [(j + 1) % L for j in range(y1, y2)]

    # Make x1 the top coordinate.
    if x1 > x2:
        x1, x2 = x2, x1

    if abs(x2 - x1) <= L - abs(x2 - x1):
        lx = abs(x2 - x1)
        x_span = [j + 1 for j in range(x1, x2)]
    else:
        lx = L - abs(x2 - x1)
        x_span = [j for j in range(x1 + 1)] + [j for j in range(x2 + 1, L)]

    # Account for the multiple paths the strings can take.
    if lx < ly:
        if abs(x2 - x1) <= L - abs(x2 - x1):
            for j in range(int((ly - lx) / 2)):
                x_span.extend(((x1 - j) % L, (x2 + j + 1) % L))
        else:
            for j in range(int((ly - lx) / 2)):
                x_span.extend(((x1 + j + 1) % L, (x2 - j) % L))
    elif lx > ly:
        for j in range(int((lx - ly) / 2)):
            y_span.extend(((y1 - j) % L, (y2 + j + 1) % L))

    y_span, x_span = list(set(y_span)), list(set(x_span))

    if ret_span:
        return (A * max([lx, ly]), y_span, x_span)
    return A * max([lx, ly])


def weight_of_edge_conditional_lattice(L, x1, y1, x2, y2, A, B1, B2,
                                       around=False, ret_span=False):
    """
    Calculates the weight of an edge between two defects in the conditional
    lattice of a concatenated glued lattice object.

    :param L: Size of LxL toric code lattice used in the concatenation.
    :type L: int
    :param x1: Row coordinate of first defect.
    :type x1: int
    :param y1: Column coordinate of first defect.
    :type y1: int
    :param x2: Row coordinate of second defect.
    :type x2: int
    :param y2: Column coordinate of second defect.
    :type y2: int
    :param B1: Weight of a horizontal or vertical edge.
    :type B1: float
    :param B2: Weight of a diagonal edge.
    :type B2: float
    :param around: Whether we should go around the lattice to match the two
    defects or not.
    :type around: bool
    :param ret_span: Whether to return the y and x coordinates spanned by the
    path between defects or not.
    :type ret_span: bool
    :return: The weight between defects and possibly the span of the path
    between them.
    :rtype: float or tuple
    """
    # If we jump between plaquettes of different color, subtract 1.
    subtract_one = False
    if (x1 + y1) % 2 != (x2 + y2) % 2:
        subtract_one = True

    # Make y1 the left coordinate.
    if y1 > y2:
        y1, y2 = y2, y1
    ly = abs(y2 - y1)
    if around:
        ly = L - abs(y2 - y1)  # ***
    y_span = [(j + 1) % L for j in range(y1, y2)]

    # Make x1 the top coordinate.
    if x1 > x2:
        x1, x2 = x2, x1

    if abs(x2 - x1) <= L - abs(x2 - x1):
        lx = abs(x2 - x1)
        x_span = [j + 1 for j in range(x1, x2)]
    else:
        lx = L - abs(x2 - x1)
        x_span = [j for j in range(x1 + 1)] + [j for j in range(x2 + 1, L)]

    if 2 * B1 <= B2:
        # Path will consist of vertical and horizontal edges only.
        wt_to_return = B1 * (ly + lx)
    else:
        # Path will consist of diagonal + vertical or diagonal + horizontal.
        diagonal_steps = min([lx, ly])
        hrz_or_vert_steps = max([lx, ly]) - diagonal_steps
        wt_to_return = B1 * hrz_or_vert_steps + B2 * diagonal_steps

    if around:
        wt_to_return += L * A  # *** To account for going on lattice1.
        y_span = [j for j in range(L)]

    if subtract_one:
        wt_to_return -= 1

    if ret_span:
        return (wt_to_return, y_span, x_span)
    return wt_to_return


def weight_of_edge(L, v, u, path_weights, cor_defc=None, is_cut=False,
                   y_cut=None, ret_interm_plaq=False, ret_span=False):
    """
    Calculates the weight of an edge between two defects located anywhere in
    a concatenated lattice object formed from two toric code lattices.

    Notes:

    * Defects take the form 'x,y,lat' where x and y are the row and column
      coordinates of the defect, respectively, and lat is lattice where it
      resides. 1 for main lattice, 2 for conditional lattice, 3 for dummy
      defect (in conditional lattice always).

    :param L: Size of LxL periodic square lattice.
    :type L: int
    :param v: First defect. In the form 'x1,y1,lat1'.
    :type v: str
    :param u: Second defect. In the form 'x2,y2,lat2'.
    :type u: str
    :param path_weights: The weights of edges in the main and conditional
    panels. In the form path_weights=(A, B1, B2) where A, B1 and B2 are the
    weights of an edge in the main panel, low weight horizontal or vertical
    edge in the conditional panel and high weight edge in the conditional
    panel, respectively.
    :type path_weights: tuple of float
    :param correc_defc: Dictionary of the 'defects to correct' associated
    with every pair of defects in the matching. Keys are the matched defects,
    (u, v), values are the associated nodes to correct, see Notes in
    func:`get_correc_defc_of_nds_rmvd`.
    :type correc_defc: dict
    :param is_cut: Whether the lattice is cut or not.
    :type is_cut: bool
    :param y_cut: The vertical coordinate of the cut, if cut.
    :type y_cut: int
    :param ret_interm_plaq: Whether to return the intermediate plaquette that
    the path which crosses a boundary traverses or not. This is only
    appropriate when one defect is in lattice 1 and the other in lattice 2.
    :type ret_interm_plaq: bool
    :param ret_span: Whether to return the y and x coordinates spanned by the
    path between defects or not.
    :type ret_span: bool
    :return: The weight between the two defects. Additionaly can return the
    intermediate boundary plaquette crossed by the path between defects and
    the coordinates of the qubits spanned by the path.
    :rtype: float or tuple
    """
    A, B1, B2 = path_weights

    x1, y1, lat1 = v.split(",")
    x2, y2, lat2 = u.split(",")
    x1, y1, lat1 = int(x1), int(y1), int(lat1)
    x2, y2, lat2 = int(x2), int(y2), int(lat2)

    # Deal with the case of dummies. 3 was used to specify them, but 2 is
    # their actual lattice placement.
    if lat1 == 3:
        lat1 = 2
    if lat2 == 3:
        lat2 = 2

    if lat1 == lat2:
        if lat1 == 1:
            assert lat2 == 1
            assert (x1 + y1) % 2 == (x2 + y2) % 2
            if cor_defc is not None:
                cor_defc.update({(v, u): (v, u), (u, v): (u, v)})
            return weight_of_edge_main_lattice(L, x1, y1, x2, y2, A,
                                               ret_span=ret_span)
        elif lat1 == 2:
            assert lat2 == 2
            if cor_defc is not None:
                cor_defc.update({(v, u): (v, u), (u, v): (u, v)})
            if not is_cut:
                return weight_of_edge_conditional_lattice(L, x1, y1,
                                                          x2, y2, A, B1, B2,
                                                          ret_span=ret_span
                                                          )
            else:
                # If both defects to one side of the cut, calculate as normal.
                if (y1 < y_cut and y2 < y_cut) or \
                   (y1 >= y_cut and y2 >= y_cut):
                    return weight_of_edge_conditional_lattice(L, x1, y1,
                                                              x2, y2, A,
                                                              B1, B2,
                                                              ret_span=ret_span
                                                              )
                # Else take the correction around the lattice.
                else:
                    return weight_of_edge_conditional_lattice(L, x1, y1,
                                                              x2, y2, A, B1,
                                                              B2, around=True,
                                                              ret_span=ret_span
                                                              )
        else:
            assert True is False
    else:
        # Make x1, y1 the coordinates of defect in lattice 1.
        if lat1 == 2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1
            u, v = v, u  # v defect in lattice 1.

        # Idea: if (x2, y2) to the right of the cut, starting at (x1, y1) we
        # must pass through the "edge" boundary. if (x2, y2) to the left of
        # the cut, starting at (x1, y1) we must pass through the "middle"
        # boundary.

        # Create dict which stores the weight of the path through each middle
        # boundary plaquette, treating each plaquette as a possible
        # intermediate step. Similarly for each edge boundary plaquette.
        middle_bd_paths, edge_bd_paths = [], []
        for x_bd in range(L):
            # Bd. plaquette must have same parity as defect in main lattice.
            if (x_bd + 0) % 2 == (x1 + y1) % 2:
                if (not is_cut) or (y2 < y_cut):
                    main_md = weight_of_edge_main_lattice(L, x1, y1, x_bd,
                                                          L, A,
                                                          ret_span=True)
                    wt_main_md, y_sp_mn, x_sp_mn = main_md
                    cond_md = weight_of_edge_conditional_lattice(L, x_bd,
                                                                 0, x2, y2,
                                                                 A, B1, B2,
                                                                 ret_span=True
                                                                 )
                    wt_cond_md, y_sp_cd, x_sp_cd = cond_md

                    wt_md = wt_main_md + wt_cond_md
                    y_sp_md = list(set(y_sp_mn + y_sp_cd))
                    x_sp_md = list(set(x_sp_mn + x_sp_cd))

                    middle_bd_paths.append((x_bd, wt_md, y_sp_md, x_sp_md))

            if (x_bd + 1) % 2 == (x1 + y1) % 2:
                if (not is_cut) or (y2 >= y_cut):
                    main_ed = weight_of_edge_main_lattice(L, x1, y1, x_bd,
                                                          -1, A,
                                                          ret_span=True)
                    wt_main_ed, y_sp_mn, x_sp_mn = main_ed
                    cond_ed = weight_of_edge_conditional_lattice(L, x_bd,
                                                                 L - 1, x2,
                                                                 y2, A, B1,
                                                                 B2,
                                                                 ret_span=True
                                                                 )
                    wt_cond_ed, y_sp_cd, x_sp_cd = cond_ed

                    wt_ed = wt_main_ed + wt_cond_ed
                    y_sp_ed = list(set(y_sp_mn + y_sp_cd))
                    x_sp_ed = list(set(x_sp_mn + x_sp_cd))

                    edge_bd_paths.append((x_bd, wt_ed, y_sp_ed, x_sp_ed))

        if len(middle_bd_paths) > 0:
            min_md_path = min(middle_bd_paths, key=itemgetter(1))
            x_bd_md, min_md, y_sp_md, x_sp_md = min_md_path
        else:
            # select something sufficiently large so that any actual weight
            # will be smaller than it.
            min_md = 2 * A * L + 1
            x_bd_md = None
        if len(edge_bd_paths) > 0:
            min_ed_path = min(edge_bd_paths, key=itemgetter(1))
            x_bd_ed, min_ed, y_sp_ed, x_sp_ed = min_ed_path
        else:
            min_ed = 2 * A * L + 1
            x_bd_ed = None

        if min_md < min_ed:
            assert x_bd_md is not None
            min_wt, x_bd, bd_type = min_md, x_bd_md, 'middle'
            y_span, x_span = y_sp_md, x_sp_md
        elif min_md > min_ed:
            assert x_bd_ed is not None
            min_wt, x_bd, bd_type = min_ed, x_bd_ed, 'edge'
            y_span, x_span = y_sp_ed, x_sp_ed
        else:  # else they are equal.
            assert (x_bd_md is not None) and (x_bd_ed is not None)
            min_wt, x_bd, bd_type = min_md, x_bd_md, 'middle'
            y_span, x_span = y_sp_md, x_sp_md

        if bd_type == 'middle':
            bd_defect = str(x_bd) + ',' + str(L) + ',1'
        else:
            bd_defect = str(x_bd) + ',' + str(-1) + ',1'

        if cor_defc is not None:
            cor_defc.update({(v, u): (v, bd_defect), (u, v): (bd_defect, v)})

        if ret_interm_plaq:
            return (min_wt, x_bd, bd_type)
        if ret_span:
            return (min_wt, y_span, x_span)

        return min_wt
