import Qubit
import Plaquette


class SquareLattice1:
    def __init__(self, size):
        """
        Initialisation of LxL periodic toric code lattice.

        Notes:

        * The dimension L of the LxL lattice must be even.
        * The lattice has alternating parity/colour plaquettes.
          Chessboard-like with dark and light plaquettes.
        * X errors light up the 2 nearby dark plaquettes.
        * Y errors light up all 4 nearby plaquettes.
        * Z errors light up the 2 nearby light plaquettes.
        * Example layout of a 4x4 lattice. The boundaries are called 'edge'
        (leftmost) and 'middle' (rightmost), which becomes useful terminology
        when concatenating lattices.

      (0,0)---(0,1)---(0,2)---(0,3)---(0,0)
        |       |       |       |       |
        |       |   .   |       |   .   |
        |       |       |       |       |
      (1,0)---(1,1)---(1,2)---(1,3)---(1,0)
        |       |       |       |       |
        |   .   |       |   .   |       |
        |       |       |       |       |
      (2,0)---(2,1)---(2,2)---(2,3)---(2,0)
        |       |       |       |       |
        |       |   .   |       |   .   |
        |       |       |       |       |
      (3,0)---(3,1)---(3,2)---(3,3)---(3,0)
        |       |       |       |       |
        |   .   |       |   .   |       |
        |       |       |       |       |
      (4,0)---(4,1)---(4,2)---(4,3)---(4,0)

        ^                               ^
        |                               |
    Edge boundary                   Middle boundary



        :param size: Dimension L of LxL lattice.
        :type size: int
        """
        assert size % 2 == 0
        self.size = size
        self.qubits = [[[] for _ in range(size)] for _ in range(size)]
        self.plaquettes = [[[] for _ in range(size)] for _ in range(size)]
        # Initialise qubits and parity check bits to trivial state.
        for i in range(size):
            for j in range(size):
                self.qubits[i][j] = Qubit.Qubit('I')
                self.plaquettes[i][j] = Plaquette.Plaquette(0)

    def apply_Y(self, i, j, boundary='none'):
        """
        Apply Y operator to qubit at position (i, j) in the lattice.

        Note:

        * Y error lights up NW, NE, SW and SE plaquettes.

        :param i: Row position of the qubit
        :type i: int
        :param j: Column position of the qubit
        :type j: int
        """
        L = self.size
        if boundary == 'middle':
            self.qubits[i][j].apply_Y()
            self.plaquettes[(i - 1) % L][(j - 1) % L].flip()  # NW
            self.plaquettes[i][(j - 1) % L].flip()  # SW
        elif boundary == 'edge':
            self.plaquettes[(i - 1) % L][j].flip()  # NE
            self.plaquettes[i][j].flip()  # SE
        elif boundary == 'none':
            L = self.size
            self.qubits[i][j].apply_Y()
            self.plaquettes[(i - 1) % L][(j - 1) % L].flip()  # NW
            self.plaquettes[(i - 1) % L][j].flip()  # NE
            self.plaquettes[i][(j - 1) % L].flip()  # SW
            self.plaquettes[i][j].flip()  # SE
        else:
            assert True is False

    def apply_Z(self, i, j, boundary='none'):
        """
        Apply Z operator to qubit at position (i,j) in the lattice.

        Note:

        * Z error lights up even parity neighbour plaquettes.

        :param i: Row position of the qubit
        :type i: int
        :param j: Column position of the qubit
        :type j: int
        """
        L = self.size
        if (i + j) % 2 == 1:
            if boundary == 'middle':
                self.qubits[i][j].apply_Z()
                self.plaquettes[i][(j - 1) % L].flip()  # SW
            elif boundary == 'edge':
                self.plaquettes[(i - 1) % L][j].flip()  # NE
            elif boundary == 'none':
                # NE and SW plaquettes are even parity.
                self.qubits[i][j].apply_Z()
                self.plaquettes[(i - 1) % L][j].flip()  # NE
                self.plaquettes[i][(j - 1) % L].flip()  # SW
            else:
                assert True is False
        else:
            if boundary == 'middle':
                self.qubits[i][j].apply_Z()
                self.plaquettes[(i - 1) % L][(j - 1) % L].flip()  # NW
            elif boundary == 'edge':
                self.plaquettes[i][j].flip()  # SE
            elif boundary == 'none':
                self.qubits[i][j].apply_Z()
                # NW and SE plaquettes are even parity.
                self.plaquettes[(i - 1) % L][(j - 1) % L].flip()  # NW
                self.plaquettes[i][j].flip()  # SE
            else:
                assert True is False

    def apply_X(self, i, j, boundary='none'):
        """
        Apply X operator to qubit at position (i,j) in the lattice.

        Note:

        * X error lights up odd parity neighbour plaquettes.

        :param i: Row position of the qubit
        :type i: int
        :param j: Column position of the qubit
        :type j: int
        """
        L = self.size
        if (i + j) % 2 == 1:
            if boundary == 'middle':
                self.qubits[i][j].apply_X()
                self.plaquettes[(i - 1) % L][(j - 1) % L].flip()  # NW
            elif boundary == 'edge':
                self.plaquettes[i][j].flip()  # SE
            elif boundary == 'none':
                self.qubits[i][j].apply_X()
                # NW and SE plaquettes are odd parity.
                self.plaquettes[(i - 1) % L][(j - 1) % L].flip()  # NW
                self.plaquettes[i][j].flip()  # SE
            else:
                assert True is False

        else:
            if boundary == 'middle':
                self.qubits[i][j].apply_X()
                self.plaquettes[i][(j - 1) % L].flip()  # SW
            elif boundary == 'edge':
                self.plaquettes[(i - 1) % L][j].flip()  # NE
            elif boundary == 'none':
                self.qubits[i][j].apply_X()
                # NE and SW plaquettes are odd parity.
                self.plaquettes[(i - 1) % L][j].flip()  # NE
                self.plaquettes[i][(j - 1) % L].flip()  # SW
            else:
                assert True is False

    def apply_stabiliser(self, stab):
        """
        Apply a stabiliser to the lattice provided a plquette coordinate.

        :param stab: The coordinates of the plquette on which the stabiliser
        should be applied. In the form stab=(x,y).
        :type stab: tuple of int
        """
        L = self.size
        x, y = stab
        corners = [(x, y), ((x + 1) % L, y), (x, (y + 1) % L),
                   ((x + 1) % L, (y + 1) % L)]
        if (x + y) % 2 == 0:
            for i, j in corners:
                self.apply_X(i, j)
        else:
            for i, j in corners:
                self.apply_Z(i, j)

    def cast_qubit_state_to_list(self):
        """
        Packages the qubit state of the lattice into a list of coordinates
        where the Pauli errors are located.

        :return: List of coordinates where Pauli X, Y and Z errors are located
        in the lattice.
        :rtype: List of list of tuples.
        """
        L = self.size
        X_coords, Y_coords, Z_coords = [], [], []
        for i in range(L):
            for j in range(L):
                q_state = self.qubits[i][j].state
                if q_state == 'X':
                    X_coords.append((i, j))
                elif q_state == 'Y':
                    Y_coords.append((i, j))
                elif q_state == 'Z':
                    Z_coords.append((i, j))
                else:
                    assert q_state == 'I'
        return [X_coords, Y_coords, Z_coords]

    def weight(self):
        """
        The weight of the Pauli operator present in the lattice.

        :return: The number of qubits where a non-trivial Pauli operator
        acts on.
        :rtype: int
        """
        L = self.size
        wt = 0
        for i in range(L):
            for j in range(L):
                if self.qubits[i][j].state != 'I':
                    wt += 1
        return wt

    def apply_correction_from_lst(self, lst):
        """
        Applies a correction operator to the lattice given a list of
        coordinates where Pauli X, Y and Z operators act on.

        :param lst: List of coordinates where X, Y and Z errors act on.
        :type lst: List of list of tuple.
        """
        X_corr, Y_corr, Z_corr = lst
        for i, j in X_corr:
            self.apply_X(i, j)
        for i, j in Y_corr:
            self.apply_Y(i, j)
        for i, j in Z_corr:
            self.apply_Z(i, j)

    def correct_in_lattice1(self, v, u, dec_meth='part_of_glued_lattice'):
        """
        Applies a correction to the lattice given the coordinates of the
        matched defects and whether the lattice is part of a concatenated
        glued lattice object or a standalone standard toric code lattice.

        :param v: Coordinate of first defect in the form 'x1,y1,lat1'. The
        final lat information dictates which lattice the defect is on. Either
        1 for main lattice, 2 for conditional lattice or 3 to indicate it's
        a dummy defect.
        :type v: str
        :param u: Coordinate of first defect in the form 'x2,y2,lat2'. The
        final lat information dictates which lattice the defect is on. Either
        1 for main lattice, 2 for conditional lattice or 3 to indicate it's
        a dummy defect.
        :type u: str
        :param dec_meth: Either 'part_of_glued_lattice', or 'standard'.
        :type dec_meth: str
        """
        L = self.size
        x1, y1, lat1 = v.split(",")
        x2, y2, lat2 = u.split(",")
        x1, y1, lat1 = int(x1), int(y1), int(lat1)
        x2, y2, lat2 = int(x2), int(y2), int(lat2)

        # Account for dummy defects whose lattice is flagged as 3 but they are
        # really in lattice 2 (the conditional lattice).
        if lat1 == 3:
            lat1 = 2
        if lat2 == 3:
            lat2 = 2

        if lat1 == lat2:
            if lat1 == 1:
                assert lat2 == 1
                assert (x1 + y1) % 2 == (x2 + y2) % 2
                if (x1 + y1) % 2 == 0:
                    # Correct along the white plaquettes.
                    if y1 <= y2:
                        self.correct_along_X_or_Z_symmetry(x1, y1, x2,
                                                           y2, dec_meth,
                                                           parity=0)
                    else:
                        self.correct_along_X_or_Z_symmetry(x2, y2, x1,
                                                           y1, dec_meth,
                                                           parity=0)
                else:
                    # Correct along the dark plaquettes.
                    if y1 <= y2:
                        self.correct_along_X_or_Z_symmetry(x1, y1, x2,
                                                           y2, dec_meth,
                                                           parity=1)
                    else:
                        self.correct_along_X_or_Z_symmetry(x2, y2, x1,
                                                           y1, dec_meth,
                                                           parity=1)
        else:
            # Make x1, y1 the coordinates of defect in lattice 1.
            if lat1 == 2:
                x1, x2 = x2, x1
                y1, y2 = y2, y1
            # Record parity of first defect, so we know along which plaquettes
            # to move in lattice 1.
            parity = (x1 + y1) % 2
            # Check whether you are matching across the middle boundary or the
            # edge boundaries
            # TODO: update these weights which clearly imply that A=B1=B2.
            bd_middle_correction_ly_wt = (L + y2) - y1
            bd_edge_correction_ly_wt = 2 * L - ((L + y2) - y1)
            if bd_middle_correction_ly_wt <= bd_edge_correction_ly_wt:
                ly_in_lat1 = L - y1
                # If defect in lattice 2 is above or at same level as defect
                # in lattice 1.
                if x2 <= x1:
                    if x1 - x2 <= L - (x1 - x2):
                        # Go up.
                        steps_up = min([ly_in_lat1, x1 - x2])
                        horizontal_steps = ly_in_lat1 - steps_up
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_up):
                            x, y = (x1 - j) % L, (y1 + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Account for potentially not moving vertically.
                        if steps_up == 0:
                            x, y = x % L, (y) % L
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)

                    else:
                        # Go down.
                        steps_down = min([ly_in_lat1, L - (x1 - x2)])
                        horizontal_steps = ly_in_lat1 - steps_down
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_down):
                            x, y = (x1 + j + 1) % L, (y1 + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Check that you moved vertically.
                        assert steps_down > 0
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)

                else:  # If defect in lattice 2 is below defect in lattice 1.
                    if x2 - x1 <= L - (x2 - x1):
                        # Go down.
                        steps_down = min([ly_in_lat1, x2 - x1])
                        horizontal_steps = ly_in_lat1 - steps_down
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_down):
                            x, y = (x1 + j + 1) % L, (y1 + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Account for potentially not moving vertically.
                        if steps_down == 0:
                            x, y = x % L, (y) % L
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                    else:
                        # Go up.
                        steps_up = min([ly_in_lat1, L - (x2 - x1)])
                        horizontal_steps = ly_in_lat1 - steps_up
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_up):
                            x, y = (x1 - j) % L, (y1 + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Check that you moved vertically.
                        assert steps_up > 0
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new + j + 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)

            elif bd_middle_correction_ly_wt > bd_edge_correction_ly_wt:
                ly_in_lat1 = y1 + 1
                # If defect in lattice 2 is above or at same level as defect
                # in lattice 1.
                if x2 <= x1:
                    if x1 - x2 <= L - (x1 - x2):
                        # Go up.
                        steps_up = min([ly_in_lat1, x1 - x2])
                        horizontal_steps = ly_in_lat1 - steps_up
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_up):
                            x, y = (x1 - j) % L, (y1 - j) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Account for potentially not moving vertically.
                        if steps_up == 0:
                            x, y = x % L, (y + 1) % L
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new - j - 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)

                    else:
                        # Go down.
                        steps_down = min([ly_in_lat1, L - (x1 - x2)])
                        horizontal_steps = ly_in_lat1 - steps_down
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_down):
                            x, y = (x1 + j + 1) % L, (y1 - j) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Check that you moved vertically.
                        assert steps_down > 0
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new - j - 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)

                else:  # If defect in lattice 2 is below defect in lattice 1.
                    if x2 - x1 <= L - (x2 - x1):
                        # Go down.
                        steps_down = min([ly_in_lat1, x2 - x1])
                        horizontal_steps = ly_in_lat1 - steps_down
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_down):
                            x, y = (x1 + j + 1) % L, (y1 - j) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Account for potentially not moving vertically.
                        if steps_down == 0:
                            x, y = x % L, (y + 1) % L
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new - j - 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                    else:
                        # Go up.
                        steps_up = min([ly_in_lat1, L - (x2 - x1)])
                        horizontal_steps = ly_in_lat1 - steps_up
                        # Apply correction, starting at x1, y1 and moving
                        # until boundary is reached.
                        x, y = x1, y1
                        for j in range(steps_up):
                            x, y = (x1 - j) % L, (y1 - j) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)
                        # Make sure that you moved vertically.
                        assert steps_up > 0
                        x_new, y_new = x, y
                        for j in range(horizontal_steps):
                            x, y = x_new % L, (y_new - j - 1) % L
                            self.apply_X(x, y) if parity == 1 else \
                                self.apply_Z(x, y)

    def correct_along_X_or_Z_symmetry(self, x1, y1, x2, y2, dec_meth, parity):
        """
        Finds and applies a correction given two defects in a toric code
        lattice which can either be standalone or as part of a concatenated
        glued lattice object.

        :param x1: Row coordinate of first defect.
        :type x1: int
        :param y1: Column coordinate of first defect.
        :type y1: int
        :param x2: Row coordinate of second defect.
        :type x2: int
        :param y2: Column coordinate of second defect.
        :param dec_meth: Either 'part_of_glued_lattice', or 'standard'.
        :type dec_meth: str
        :param parity: Determines which symmetry we move along. 0 for Z
        symmetry. 1 for X symmetry.
        :type parity: int
        """
        # Order is important! Start at (x1, y1) and match inside the lattice
        # (horizontally) to get to (x2, y2).
        L = self.size

        # Check if we are crossing boundaries.
        bd = None
        if y1 == -1:
            bd = 'edge'
        elif y2 == L:
            bd = 'middle'

        # Check that second defect is to the right of first defect.
        assert y2 >= y1
        if parity == 1:
            # Check that points lie on X symmetry.
            assert (x1 + y1) % 2 == 1
            assert (x2 + y2) % 2 == 1
        else:
            # Check that points lie on Z symmetry.
            assert (x1 + y1) % 2 == 0
            assert (x2 + y2) % 2 == 0
        # print(x1, y1, x2, y2, bd)
        if (y2 - y1 <= L - (y2 - y1)) or dec_meth == 'part_of_glued_lattice':
            # Correct inside the lattice.
            ly = y2 - y1  # Horizontal distance between defects.
            # If 2nd defect is above or at same level as 1st defect.
            if x2 <= x1:
                if x1 - x2 <= L - (x1 - x2):
                    lx = x1 - x2
                    # Go up.
                    steps_up = min([ly, lx])
                    # Apply correction, start at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    # Keep track of where you first applied a gate.
                    x_ini, y_ini = None, None
                    for j in range(steps_up):
                        x, y = (x1 - j) % L, (y1 + j + 1) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                        if x_ini is None:
                            x_ini, y_ini = x, y
                    if lx < ly and ((bd is None) or (bd == 'edge') or ((bd == 'middle') and ((y2 - y) % L > 0))):
                        # Transform back from qubit coord to plaquette coord if moved.
                        if steps_up > 0:
                            x, y = (x - 1) % L, y % L
                        x, y, x0, y0 = self.correct_horizontally(x, y, x2, y2, parity)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    elif lx >= ly and min([abs(x2 - x1), L - abs(x2 - x1)]) > 0:
                        # Transform back from qubit coord to plaquette coord if moved.
                        if steps_up > 0:
                            x, y = (x - 1) % L, y % L
                        x, y, x0, y0 = self.correct_vertically(x, y, x2, y2, parity, up=True)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    x_md, y_md = x, y
                    x_ed, y_ed = x_ini, y_ini
                else:
                    lx = L - (x1 - x2)
                    # Go down.
                    steps_down = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    x_ini, y_ini = None, None  # Keep track of where you first applied a gate.
                    for j in range(steps_down):
                        x, y = (x1 + j + 1) % L, (y1 + j + 1) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                        if x_ini is None:
                            x_ini, y_ini = x, y
                    if lx < ly and ((bd is None) or (bd == 'edge') or ((bd == 'middle') and ((y2 - y) % L > 0))):
                        # Transform back from qubit coord to plaquette coord if moved.
                        if steps_down > 0:
                            x, y = x % L, y % L
                        x, y, x0, y0 = self.correct_horizontally(x, y, x2, y2, parity)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    elif lx >= ly and min([abs(x2 - x1), L - abs(x2 - x1)]) > 0:
                        # Transform back from qubit coord to plaquette coord if moved.
                        if steps_down > 0:
                            x, y = x % L, y % L
                        x, y, x0, y0 = self.correct_vertically(x, y, x2, y2, parity, up=False)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    x_md, y_md = x, y
                    x_ed, y_ed = x_ini, y_ini

            else:  # If 2nd defect is below 1st defect.
                if x2 - x1 <= L - (x2 - x1):
                    lx = x2 - x1
                    # Go down.
                    steps_down = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    x_ini, y_ini = None, None  # Keep track of where you first applied a gate.
                    for j in range(steps_down):
                        x, y = (x1 + j + 1) % L, (y1 + j + 1) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                        if x_ini is None:
                            x_ini, y_ini = x, y
                    if lx < ly and ((bd is None) or (bd == 'edge') or ((bd == 'middle') and ((y2 - y) % L > 0))):
                        # Transform back from qubit coord to plaquette coord if moved
                        if steps_down > 0:
                            x, y = x % L, y % L
                        x, y, x0, y0 = self.correct_horizontally(x, y, x2, y2, parity)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    elif lx >= ly and min([abs(x2 - x1), L - abs(x2 - x1)]) > 0:
                        # Transform back from qubit coord to plaquette coord if moved
                        if steps_down > 0:
                            x, y = x % L, y % L
                        x, y, x0, y0 = self.correct_vertically(x, y, x2, y2, parity, up=False)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    x_md, y_md = x, y
                    x_ed, y_ed = x_ini, y_ini

                else:
                    lx = L - (x2 - x1)
                    # Go up.
                    steps_up = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    x_ini, y_ini = None, None  # Keep track of where you first applied a gate.
                    for j in range(steps_up):
                        x, y = (x1 - j) % L, (y1 + j + 1) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                        if x_ini is None:
                            x_ini, y_ini = x, y
                    if lx < ly and ((bd is None) or (bd == 'edge') or ((bd == 'middle') and ((y2 - y) % L > 0))):
                        # Transform back from qubit coord to plaquette coord if moved.
                        if steps_up > 0:
                            x, y = (x - 1) % L, y % L
                        x, y, x0, y0 = self.correct_horizontally(x, y, x2, y2, parity)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    elif lx >= ly and min([abs(x2 - x1), L - abs(x2 - x1)]) > 0:
                        # Transform back from qubit coord to plaquette coord if moved.
                        if steps_up > 0:
                            x, y = (x - 1) % L, y % L
                        x, y, x0, y0 = self.correct_vertically(x, y, x2, y2, parity, up=True)
                        if x_ini is None:
                            x_ini, y_ini = x0, y0
                    x_md, y_md = x, y
                    x_ed, y_ed = x_ini, y_ini

            if bd == 'middle':
                # Undo last correction and apply appropriate boundary
                # operator.
                self.apply_X(x_md, y_md) if parity == 1 else self.apply_Z(x_md, y_md)
                self.apply_X(x_md, y_md, bd) if parity == 1 else \
                    self.apply_Z(x_md, y_md, bd)
            elif bd == 'edge':
                # Undo first correction and apply appropriate boundary
                # operator.
                assert x_ed is not None and y_ed is not None
                self.apply_X(x_ed, y_ed) if parity == 1 else self.apply_Z(x_ed, y_ed)
                self.apply_X(x_ed, y_ed, bd) if parity == 1 else \
                    self.apply_Z(x_ed, y_ed, bd)
            else:
                assert bd is None

        elif (y2 - y1 > L - (y2 - y1)):
            assert dec_meth == 'standard'
            # Correct around the lattice.
            # TODO: Fill this method.
            ly = L - (y2 - y1)  # Horizontal distance between defects.
            # If 2nd defect is above or at same level as 1st defect.
            if x2 <= x1:
                if x1 - x2 <= L - (x1 - x2):
                    lx = x1 - x2
                    # Go up.
                    steps_up = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    for j in range(steps_up):
                        x, y = (x1 - j) % L, (y1 - j) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                    # Transform back from qubit coord to plaquette coord if moved.
                    if steps_up > 0:
                        x, y = (x - 1) % L, (y - 1) % L
                    if lx < ly:
                        self.correct_horizontally_to_left(x, y, x2, y2, parity)
                    else:
                        self.correct_vertically(x, y, x2, y2, parity, up=True)

                else:
                    lx = L - (x1 - x2)
                    # Go down.
                    steps_down = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    for j in range(steps_down):
                        x, y = (x1 + j + 1) % L, (y1 - j) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                    # Transform back from qubit coord to plaquette coord if moved.
                    if steps_down > 0:
                        x, y = x % L, (y - 1) % L
                    if lx < ly:
                        self.correct_horizontally_to_left(x, y, x2, y2, parity)
                    else:
                        self.correct_vertically(x, y, x2, y2, parity, up=False)

            else:  # If 2nd defect is below 1st defect.
                if x2 - x1 <= L - (x2 - x1):
                    lx = x2 - x1
                    # Go down.
                    steps_down = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    for j in range(steps_down):
                        x, y = (x1 + j + 1) % L, (y1 - j) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                    # Transform back from qubit coord to plaquette coord if moved
                    if steps_down > 0:
                        x, y = x % L, (y - 1) % L
                    if lx < ly:
                        self.correct_horizontally_to_left(x, y, x2, y2, parity)
                    else:
                        self.correct_vertically(x, y, x2, y2, parity, up=False)

                else:
                    lx = L - (x2 - x1)
                    # Go up.
                    steps_up = min([ly, lx])
                    # Apply correction, starting at x1, y1 and ending at x2, y2.
                    x, y = x1, y1
                    for j in range(steps_up):
                        x, y = (x1 - j) % L, (y1 - j) % L
                        self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                    # Transform back from qubit coord to plaquette coord if moved.
                    if steps_up > 0:
                        x, y = (x - 1) % L, (y - 1) % L
                    if lx < ly:
                        self.correct_horizontally_to_left(x, y, x2, y2, parity)
                    else:
                        self.correct_vertically(x, y, x2, y2, parity, up=True)

    def correct_horizontally(self, x1, y1, x2, y2, parity):
        """
        Finds and applies a correction between two defects which are on the
        same row.

        :param x1: Row coordinate of first defect.
        :type x1: int
        :param y1: Column coordinate of first defect.
        :type y1: int
        :param x2: Row coordinate of second defect.
        :type x2: int
        :param y2: Column coordinate of second defect.
        :param parity: Determines which symmetry we move along. 0 for Z
        symmetry. 1 for X symmetry.
        :type parity: int
        :return: The coordinates where the first gate is applied, (x0,y0), and
        the coordinates where the last gate is applied, (x,y).
        :rtype: tuple
        """
        L = self.size
        x0, y0 = None, None
        for j in range(y2 - y1):
            x, y = x1 % L, (y1 + 1 + j) % L
            self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
            if x0 is None:
                x0, y0 = x, y
        return (x, y, x0, y0)

    def correct_horizontally_to_left(self, x1, y1, x2, y2, parity):
        """
        Finds and applies a correction between two defects which are on the
        same row. The correction is applied starting at the first defect and
        moving left to the second defect, potentially around the lattice.

        :param x1: Row coordinate of first defect.
        :type x1: int
        :param y1: Column coordinate of first defect.
        :type y1: int
        :param x2: Row coordinate of second defect.
        :type x2: int
        :param y2: Column coordinate of second defect.
        :param parity: Determines which symmetry we move along. 0 for Z
        symmetry. 1 for X symmetry.
        :type parity: int
        """
        L = self.size
        if y1 >= y2:
            for j in range((y1 - y2) % L):
                x, y = x1 % L, (y1 - j) % L
                self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
        else:
            for j in range(L - (y2 - y1)):
                x, y = x1 % L, (y1 - j) % L
                self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)

    def correct_vertically(self, x1, y1, x2, y2, parity, up):
        """
        Finds and applies a correction between two defects which are on the
        same column.

        :param x1: Row coordinate of first defect.
        :type x1: int
        :param y1: Column coordinate of first defect.
        :type y1: int
        :param x2: Row coordinate of second defect.
        :type x2: int
        :param y2: Column coordinate of second defect.
        :param parity: Determines which symmetry we move along. 0 for Z
        symmetry. 1 for X symmetry.
        :type parity: int
        :param up: Whether to move up starting at the first defect or not
        (i.e. move down).
        :type up: bool
        :return: The coordinates where the first gate is applied, (x0,y0), and
        the coordinates where the last gate is applied, (x,y).
        :rtype: tuple
        """
        L = self.size
        x, y = x1, y1
        x0, y0 = None, None
        if up:
            for j in range(min([abs(x2 - x1), L - abs(x2 - x1)])):
                x, y = (x1 - j) % L, y1 % L
                self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                if x0 is None:
                    x0, y0 = x, y
        else:
            for j in range(min([abs(x2 - x1), L - abs(x2 - x1)])):
                x, y = (x1 + j + 1) % L, y1 % L
                self.apply_X(x, y) if parity == 1 else self.apply_Z(x, y)
                if x0 is None:
                    x0, y0 = x, y
        return (x, y, x0, y0)

    def is_in_trivial_state_X1(self):
        """
        Checks for type 1 logical operators. These are rows of X operators
        (XXXX) or any muplication of these by stabilisers. Ignore Z errors.

        :return: Whether it is in trivial state with respect to the logical
        operators or not.
        :rtype: Bool
        """
        L = self.size
        # Check parity along columns
        # Note we exclude boundary qubits since the applied correction is true
        # up to the boundaries and there are remnant defects there.
        for j in range(1, L - 1):
            total = 0
            for i in range(L):
                if self.qubits[i][j].state == 'Y' or \
                     self.qubits[i][j].state == 'X':
                    total += 1
            if total % 2 == 1:
                return False
        return True

    def is_in_code_space(self):
        """
        Checks whether the state of the lattice is in the code space and
        therefore suitable for making an inference about its logical error
        state.

        :return: Whether it is in the code space or not.
        :rtype: Bool
        """
        L = self.size
        for i in range(L):
            for j in range(L):
                if self.plaquettes[i][j].state == 1:
                    return False
        return True

    def are_defects_at_boundary(self):
        """
        Checks whether the defects on the lattice are at the boundary
        plaquettes, in the horizontal direction.

        Notes:

        * This method should be used to check that the remnant defects all lie
          at the boundary after a correction is applied when the lattice is
          part of a concatenated glued lattice object.

        :return: Whether the defects in the lattice or lie at the boundary
        plaquettes or not.
        :rtype: Bool
        """
        L = self.size
        for i in range(L):
            for j in range(L):
                if (j != 0) and (j != L - 1) and (self.plaquettes[i][j].state == 1):
                    return False
        return True

    def print_plaquettes(self):
        """
        Prints the plaquette state of the lattice, i.e. where the defects lie.
        """
        for i in range(self.size):
            for j in range(self.size):
                print(self.plaquettes[i][j].state, end="")
            print()
        print()

    def __repr__(self):
        strn = []
        for i in range(self.size):
            for j in range(self.size):
                strn.append(self.qubits[i][j].state)
            strn.append('\n')
        return "".join(strn)
