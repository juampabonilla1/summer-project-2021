class Qubit:
    """
    Defines a physical qubit and its state under Pauli noise

    Notes:

    * This is used as a helper class to the 'Lattice' class.

    Use cases:

    * Get qubit error state: :meth:`print_qubit`.
    * Apply Pauli gate to qubit: :meth:`apply_Z`, :meth:`apply_X`,
    :meth:`apply_Y`.
    """

    def __init__(self, state):
        """
        Initialisation of qubit.

        :param state: State of qubit under Pauli noise. Either 'X', 'Y', 'Z'
        or 'I' (trivial state).
        :type state: string
        """
        self.state = state

    def apply_Z(self):
        """
        Apply Pauli Z operator to qubit.
        """
        if self.state == 'I':
            self.state = 'Z'
        elif self.state == 'X':
            self.state = 'Y'
        elif self.state == 'Y':
            self.state = 'X'
        elif self.state == 'Z':
            self.state = 'I'

    def apply_X(self):
        """
        Apply Pauli X operator to qubit.
        """
        if self.state == 'I':
            self.state = 'X'
        elif self.state == 'X':
            self.state = 'I'
        elif self.state == 'Y':
            self.state = 'Z'
        elif self.state == 'Z':
            self.state = 'Y'

    def apply_Y(self):
        """
        Apply Pauli Y operator to qubit.
        """
        if self.state == 'I':
            self.state = 'Y'
        elif self.state == 'X':
            self.state = 'Z'
        elif self.state == 'Y':
            self.state = 'I'
        elif self.state == 'Z':
            self.state = 'X'

    def __repr__(self):
        return ''.join(self.state)
