def is_coord_in_tuple(defect, tuple):
    """
    Checks whether a certain defect coordinate is in a tuple of defects.
    Checks for coordinate equality, disregarding the lattice the defect is in.

    :param defect: Defect in the form 'x,y,lat', where x,y are the desired
    target coordinates.
    :type defect: str
    :param tuple: Tuple of defects in the form (u, v).
    :type tuple: tuple
    :return: Whether defect coordinates are in the tuple or not.
    :rtype: bool
    """
    defect1, defect2 = tuple
    x, y, lat = defect.split(",")
    x1, y1, lat1 = defect1.split(",")
    x2, y2, lat2 = defect2.split(",")
    standard = ['1', '2']
    if lat in standard:
        cond1 = ((x, y) == (x1, y1) and lat1 in standard)
        cond2 = ((x, y) == (x2, y2) and lat2 in standard)
        return cond1 or cond2
    else:  # Else deal with dummy equality.
        assert lat == '3'  # flag to identify dummy defects.
        cond1 = (x, y, lat) == (x1, y1, lat1)
        cond2 = (x, y, lat) == (x2, y2, lat2)
        return cond1 or cond2


def are_coords_equal(defect1, defect2):
    """
    Checks for equality of defect coordinates, disregarding the lattice they
    are on.

    :param defect1: The first defect, in the form 'x1,y1,lat1'.
    :type defect1: str
    :param defect2: The second defect, in the form 'x2,y2,lat2'.
    :type defect2: str
    :return: Whether the two defects have the same coordinates or not.
    :rtype: bool
    """
    x1, y1, lat1 = defect1.split(",")
    x2, y2, lat2 = defect2.split(",")
    if lat1 in ['1', '2'] and lat2 in ['1', '2']:
        return (x1, y1) == (x2, y2)
    else:  # Else deal with dummy equality. lat=3 used to flag dummy defects.
        assert (lat1 == '3') or (lat2 == '3')
        return (x1, y1, lat1) == (x2, y2, lat2)
