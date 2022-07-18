from numpy.random import binomial as bino
from qiskit import Aer
from qiskit import transpile

def sample(means, arm_no):
	return bino(1, means[arm_no])

def qc_to_unitary(qc):
	backend = Aer.get_backend('unitary_simulator')
	transpiled = transpile(qc, backend = backend)
	return backend.run(transpiled).result().get_unitary()

def get_statevector(qc):
	backend = Aer.get_backend('statevector_simulator')
	transpiled = transpile(qc, backend = backend)
	return backend.run(transpiled).result().get_statevector()