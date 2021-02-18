class Plaquette:
    """
    Defines a parity check on a lattice plaquette. For the XZZX surface code,
    physical qubits are represented as vertices and parity check bits can be
    considered to lie in the plaquettes.

    Notes:

    * This is used as a helper class to the 'Lattice' class.

    Use cases:

    * Get plaquette state: :meth:`print_qubit`.
    * Change the plaquette state: :meth:`flip`.
    """

    def __init__(self, state):
        """
        Initialisation of plaquette check bit.

        :param state: State of parity check plaquette. Either 1 or 0 (trivial
        state).
        :type state: int
        """
        self.state = state

    def flip(self):
        """
        Chage the state of the read-out bit.
        """
        self.state = self.state ^ 1

    def __repr__(self):
        return ''.join(str(self.state))
