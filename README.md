# summer-project-2021

Python software to model quantum error correction by concatenation of toric code lattices.

## Installation

```bash
$ git clone https://github.com/juampabonilla1/summer-project-2021.git
```

## Usage

The program [src/simulation_starter.py](src/simulation_starter.py) runs N simulations of error input, syndrome measurement, recovery operator creation and error correction on the square periodic toric code lattice. The noise model is one that is biased toward dephasing, see [arXiv:1708.08474](https://arxiv.org/abs/1708.08474), [arXiv:1907.02554](https://arxiv.org/abs/1907.02554), [arXiv:1812.08186](https://arxiv.org/abs/1812.08186), [arXiv:2009.07851](https://arxiv.org/abs/2009.07851). When bias tends to infinity, noise tends to purely dephasing. When bias equals 0.5, noise equals to the depolarising channel.

  * Example: Run 100 rounds of Monte Carlo simulations on a 6x6 toric
    code lattice with a noise model with physical error probability of
    0.15 on every qubit. The bias of the noise toward dephasing is 0.5,
    so that this is the depolarising channel.

```bash
python src/simulation_starter.py p=0.15 L=6 N=100 bias=0.5
...
[{'n_fail: 19', 'measurement_error_probability: 0', 'decoder: Rotated square XZ MWPM', 'code: Rotated square XZ 6x6', 'n_run: 100', 'bias: 0.5', 'error_model: Biased noise toward dephasing', 'logical_failure_rate: 0.19', 'wall_time: 3.406278610229492', 'error_probability: 0.15'}]
```

### Matching libraries

In order to perform minimum weight perfect matching (the bottleneck task) using an optimised program, install [qecsim](https://qecsim.github.io/).

After a successful installation, comment out the following lines in [src/graph_functions.py](src/graph_functions.py):

```python
import networkx as nx
import matplotlib.pyplot as plt
import minimum_weight_perfect_matching as mwpm
```

Uncomment the following lines in [src/graph_functions.py](src/graph_functions.py):

```python
import qecsim.graphtools as gt
```

The program can now be run as normal and a fast implementation of MWPM will be used.

## Acknowledgments

This project makes use of the tools provided by [qecsim](https://qecsim.github.io/):

D. K. Tuckett, Tailoring surface codes: Improvements in quantum error correction with biased noise, [Ph.D. thesis](https://hdl.handle.net/2123/22132), University of Sydney (2020), (qecsim: https://github.com/qecsim/qecsim)

## License
[MIT](https://choosealicense.com/licenses/mit/)