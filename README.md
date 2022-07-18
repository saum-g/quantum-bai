## A High-level Simulation of Quantum Algorithms for Best Arm Identification in Quantum Multi-Armed Bandits

This is a high-level simulation of the algorithm presented by [Wang et al](https://arxiv.org/abs/2007.07049v2) to calculate the number of calls made to the QAB oracle by the algorithm. 

In order to call the algorithm for a specific bandit instance, call the `BestArm` function in the [abstract_quantum_simulation.py](abstract_quantum_simulation.py) file. The files [get_results_H.py](get_results_H.py) and [get_results_delta.py](get_results_delta.py) are convenience files to derive statistics over various 2-armed bandit instances.
