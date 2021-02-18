import random

import SquareLattice1 as sl


class GluedLattices:
    def __init__(self, size):
        """
        Defines a concatenated Lx(2L) glued lattice from two standard LxL
        toric code lattices. Also creates a separate standard LxL toric code
        lattice.

        Notes:

        * The separate LxL toric code lattice is created to investigate the
        effect of an equal error model on the glued lattice object compared
        to on the standard toric code lattice.

        :param size: Dimension L of LxL toric code lattices used for
        concatenation.
        :type size: int
        """
        self.size = size
        self.lattice1 = sl.SquareLattice1(size)
        self.lattice2 = sl.SquareLattice1(size)
        self.lattice_standard = sl.SquareLattice1(size)

    def add_high_rate_Z_noise(self, p, n):
        """
        Add biased dephasing (Z) Pauli noise. The same noise is added to each
        lattice.

        :param p: Probability of applying an error to any one qubit.
        :type p: float
        :param n: Bias coefficient
        :type n: float
        """
        L = self.size
        # Calculate probabilities of each of a Z, X and Y error.
        pz = (p * n) / (n + 1)
        px = (p) / (2 * (n + 1))
        py = (p) / (2 * (n + 1))
        # Normalise to the probability of an error, p
        normaliser = 1 / p
        pz_n = pz * normaliser
        py_n = py * normaliser
        px_n = px * normaliser
        for i in range(L):
            for j in range(L):
                # Check if error will be applied to the qubit at index (i,j).
                if random.random() <= p:
                    # Check what type of error will be applied.
                    r = random.random()
                    if r <= pz_n:
                        self.lattice1.apply_Z(i, j)
                        self.lattice2.apply_Z(i, j)
                        self.lattice_standard.apply_Z(i, j)
                    elif pz_n < r <= (pz_n + px_n):
                        self.lattice1.apply_X(i, j)
                        self.lattice2.apply_X(i, j)
                        self.lattice_standard.apply_X(i, j)
                    elif (px_n + pz_n) < r <= (pz_n + px_n + py_n):
                        self.lattice1.apply_Y(i, j)
                        self.lattice2.apply_Y(i, j)
                        self.lattice_standard.apply_Y(i, j)

    def __repr__(self):
        L = self.size
        strn = []
        for i in range(L):
            for j in range(L):
                strn.append(self.lattice1.qubits[i][j].state)
            strn.append('\n')
        strn.append('\n')
        for i in range(L):
            for j in range(L):
                strn.append(self.lattice2.qubits[i][j].state)
            strn.append('\n')
        strn.append('\n')
        for i in range(L):
            for j in range(L):
                strn.append(self.lattice_standard.qubits[i][j].state)
            strn.append('\n')
        return "".join(strn)
