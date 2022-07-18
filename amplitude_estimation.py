from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
from qiskit import Aer
from qiskit.tools.visualization import circuit_drawer
from qiskit.quantum_info import Operator
from qiskit.circuit.library import QFT, IntegerComparator
from qiskit.circuit.library.standard_gates import XGate
from util import *
import math
import numpy as np

def ith_arm_oracle(i, means):
	# return unitary O_p_i

	pi = means[i]
	unit = [[math.sqrt(1-pi), math.sqrt(pi)], 
	[math.sqrt(pi), -math.sqrt(1-pi)]]

	return Operator(unit)

def get_Q(prep_unit, no_state_qubits):
	# return unitary gate corresponding to Q = -AS_0A^-1S_X
	# assume last qubit in state is the flag

	q = QuantumRegister(no_state_qubits)
	qc = QuantumCircuit(q)

	# if last is 1 then put -1
	qc.z(q[-1])

	# unitary matrix => inverse is adjoint
	prep_inv = prep_unit.adjoint()
	qc.append(prep_inv, qargs = q)

	# S0 changes sign of state 0
	s0_unit = np.eye(2**no_state_qubits)
	s0_unit[0,0] = -1
	s0_op = Operator(s0_unit)
	qc.append(s0_op, qargs = q)

	qc.append(prep_unit, qargs = q)

	Q = qc_to_unitary(qc)

	# multiply by -1
	return Operator(Q.data*-1).to_instruction()

def get_qft(no_qbits, inverse = False):
	return qc_to_unitary(QFT(no_qbits, inverse = inverse))


def amplitude_est(prep_unit, no_state_qubits, no_prec_qubits, sim = 'unitary_simulator', num_shots = 1):
	"""
	Perform amplitude estimation on state prepared by prep_unit
	Input:
	no_state_qubits: number of qubits representing the state
	assume last qubit in state is the flag
	1/2^no_prec_qubits precision on value
	num_shots: no of samples to derive if aer_simulator used
	Output:
	no of calls to prep_unit AND
	post-processing result of amplitude estimation depending on the simulator used
	"""

	if 'simulator_' in sim:
		backend = provider.get_backend(sim)

	else:
		backend = Aer.get_backend(sim)
	# backend = Aer.get_backend('aer_simulator')
	# backend = Aer.get_backend('statevector_simulator')
	# backend = Aer.get_backend('unitary_simulator')

	s = QuantumRegister(no_state_qubits)
	p = QuantumRegister(no_prec_qubits)
	if sim == 'aer_simulator':
		res = ClassicalRegister(no_prec_qubits)
		qc = QuantumCircuit(s, p, res)
	else:
		qc = QuantumCircuit(s, p)

	qc.append(prep_unit, qargs = s)

	# return transpile(qc, backend = backend)

	# apply F_M
	F_M = get_qft(no_prec_qubits)
	# qc.barrier()
	qc.append(F_M, qargs = p)
	# qc.barrier()

	# return transpile(qc, backend = backend)

	Q = get_Q(prep_unit, no_state_qubits)

	# Gamma_M(Q)
	for i in range(no_prec_qubits):
		exponent = 2**(i)
		label = "Q^("+str(exponent)+")"
		Q_exp_j = Q.power(exponent).control(1, label = label)

		# if ith qubit is 1, apply Q^(2^(i))
		qc.append(Q_exp_j, qargs = [p[i]]+list(s))

	# return transpile(qc, backend = backend)

	# F_M-1
	F_M_inv = get_qft(no_prec_qubits, inverse = True)
	qc.append(F_M_inv, qargs = p)


	if sim == 'aer_simulator':
		qc.measure(p, res)

	transpiled = transpile(qc, backend = backend)

	# print(transpiled.draw())

	no_calls_unit = 2**(no_prec_qubits + 1) - 1

	if 'statevector' in sim:
		return no_calls_unit, backend.run(transpiled).result().get_statevector()
	elif sim == 'aer_simulator':
		return no_calls_unit, backend.run(transpiled, shots = num_shots).result().get_counts()
	else:
		return no_calls_unit, qc_to_unitary(transpiled)
	# print(backend.run(transpiled).result().get_statevector())

	# res = backend.run(transpiled, shots = 1000).result().get_counts()

	# print(res)

	# for i in res.keys():
	# 	if res[i] > 0:
	# 		return math.sin(math.pi*int(i,2)/(2**no_prec_qubits))**2

	# return qc_to_unitary(transpiled)

