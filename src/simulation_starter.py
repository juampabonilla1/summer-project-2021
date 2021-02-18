"""
When running on the cluster, in module `graph_functions`, comment out:
--> 'import networkx as nx'
--> 'import matplotlib.pyplot as plt'
--> 'import minimum_weight_perfect_matching as mwpm'

and uncomment:
--> 'import qecsim.graphtools as gt'
"""

import math
import time
import sys

import GluedLattices as gl
from glued_lattice_decoder import decode_glued_lattice
from standard_toric_code_decoder import standard_decode
# import error_configurations_to_correct as ec


def test(p, bias, size, N):
    """
    Run N simulations of error input and error correction on a square lattice.

    :param p: Probability of an error on any one physical qubit.
    :type p: float
    :param bias: Bias coefficient. float greater than 0 quantifying the
    bias toward dephasing.
    :type bias: float
    :param size: Size L of LxL lattice. Must be an even number.
    :type size: int
    :param N: Number of simulations to run.
    :type N: int
    """
    # TODO: create weights and store in path_weights = (A, B1, B2)
    # where A is the weight of a step in lattice1. B1 and B2 are the weights
    # of a high rate (horizontal) and low rate (diagonal) steps in lattice2,
    # respectively.
    B = -math.log(0.5)
    A = -math.log((2 * p / 3) / (1 - p))
    path_weights = (A, B, B)

    path_weights = (1, 1, 1)

    start_time = time.time()  # Time simulations.

    assert size % 2 == 0

    # Simulate errors and correction.
    total_fail = 0
    for _ in range(N):
        # print("__________")
        lattice = gl.GluedLattices(size)
        lattice.add_high_rate_Z_noise(p, bias)

        # print("#########################")
        # print(lattice.lattice1)
        # print(lattice.lattice_standard)

        decode_glued_lattice(lattice, path_weights)
        # standard_decode(lattice, p, bias)

        # print(lattice.lattice1)
        # print(lattice.lattice_standard)
        # print("$$$$$$$$$$$$$$$$$$$$$$")

        # if not lattice.lattice_standard.is_in_code_space():
        #    assert True is False, "Standard lattice not in code space!"

        if not lattice.lattice1.are_defects_at_boundary():
            print(lattice.lattice1)
            lattice.lattice1.print_plaquettes()
            assert True is False, "Remnant defects should be at boundary"

        # standard_fail, glued_fail = False, False
        if (not lattice.lattice1.is_in_trivial_state_X1()):
            # glued_fail = True
            total_fail += 1
        # if (not lattice.lattice_standard.is_in_trivial_state_X1()):
        #    standard_fail = True
        # if glued_fail and not standard_fail:
        #    print('Check this one out!')
        #    return

    print([{"code: Rotated square XZ " + str(size) + "x" + str(size),
            "decoder: Rotated square XZ MWPM",
            "error_model: Biased noise toward dephasing",
            "bias: " + str(bias),
            "error_probability: " + str(p),
            "logical_failure_rate: " + str(total_fail / N),
            "measurement_error_probability: " + str(0),
            "n_run: " + str(N),
            "n_fail: " + str(total_fail),
            "wall_time: " + str(time.time() - start_time)}])

    # print(bias, size, p, total_fail, N, time.time() - start_time)


def start():
    try:
        assert len(sys.argv) >= 5
        all_arguments = []
        for i in range(len(sys.argv)):
            if i > 0:
                argument = sys.argv[i]
                key, value = argument.split('=')
                if key == "p":
                    all_arguments.append("p")
                    p = float(value)
                elif key == "bias":
                    all_arguments.append("bias")
                    if value == 'inf':
                        bias = 'inf'
                    else:
                        bias = float(value)
                elif key == "L":
                    all_arguments.append("L")
                    L = int(value)
                    assert L % 2 == 0
                elif key == "N":
                    all_arguments.append("N")
                    N = int(value)
        assert all(item in all_arguments for item in ["p", "bias", "L", "N"])

    except Exception:
        print("#######################")
        print()
        print("Incorrect input syntax.")
        print(".......................")
        print("Run this program as: ")
        print()
        print("python square_periodic_lattice_decoder.py"
              + " p=<physical_error_rate>"
              + " bias=<bias_coefficient> L=<size_of_square_lattice>"
              + " N=<tot_number_of_simulations>")
        print()
        print("<physical_error_rate>: float between 0 and 1 (inclusive)")
        print("<bias_coefficient>: float greater than 0. Can also be str"
              + " 'inf' to use a noise model infinitely biased toward"
              + " dephasing")
        print("<size_of_square_lattice>: int positive and even.")
        print("<tot_number_of_simulations>: int greater than zero.")
        print(".......................")
        print("Example: ")
        print("python square_periodic_lattice_decoder.py"
              + " p=0.25 bias=10 L=24 N=10")
        print()
        print("#######################")

        sys.exit()

    test(p, bias, L, N)


if __name__ == "__main__":
    start()

# test(0.155, 0.5, 4, 1)
# test(0.155, 0.5, 6, 1000)
# test(0.155, 0.5, 8, 1000)
