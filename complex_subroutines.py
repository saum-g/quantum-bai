import math
import statistics
import random
from amplitude_estimation import *

def mult_bound_estimate(prob_1, c, k, p = 0.01):
	"""
	Outputs:
	number of calls to oracle preparing state with prob prob_1 AND
	estimate of prob_1 upto c multiplicative error if prob_1>=p, 0 if prob_1=0 with probability >=1-1/2^k
	"""
	M_max = (8 * math.pi)/(c * math.sqrt((1 - c) * p))
	# print("M_max =", M_max)

	no_prec_qubits = 1
	M = 2**no_prec_qubits

	total_calls_to_oracle = 0

	while True:
		# print("\n######")
		if M > M_max:
			# print("returning after reaching very high M")
			return total_calls_to_oracle, 0

		reps = math.ceil(((2 + k + math.log(math.log(M_max, 2), 2)))*math.log(2) / (2 * ((8/(math.pi**2)) - 1/2)))
		calls, ae_counts = amplitude_est(ith_arm_oracle(0, [prob_1]),
				1, no_prec_qubits, sim = 'aer_simulator', num_shots = reps)

		# expt was repeated reps times
		calls *= reps

		total_calls_to_oracle += calls

		# print("reps=", reps)
		# print("counts=", ae_counts)

		estimates = []
		for val, shots in ae_counts.items():
			est = math.sin(math.pi * int(val, 2)/(M))**2
			estimates += [est]*shots

		# print("estimates found:", estimates)

		med_est = statistics.median(estimates)
		# print("median:", med_est)

		error = 2 * math.pi * math.sqrt(med_est * (1 - med_est))/M + (math.pi ** 2)/(M**2)
		# print("error=", error)

		if  error <= c * med_est:
			return total_calls_to_oracle, med_est

		else:
			# print("error higher than threshold=", c * med_est)
			no_prec_qubits += 1
			M *= 2

def search(prob_1):
	"""
	Output:
	number of calls to oracle for computing the state till correct value is found
	"""
	l = 0
	c = random.uniform(1.001, 1.999)

	no_calls = 0

	while True:
		l += 1
		M = math.ceil(c**l)
		sample = random.choices([0, 1], [1 - prob_1, prob_1])[0]
		no_calls += 1
		if sample == 1:
			return no_calls
		j = random.randint(1, M)

		q = QuantumRegister(1)
		res = ClassicalRegister(1)
		qc = QuantumCircuit(q, res)

		unit = ith_arm_oracle(0, [prob_1])
		Q = get_Q(unit, 1)

		qc.append(unit, qargs = q)
		no_calls += 1

		qc.append(Q.power(j), qargs = q)
		no_calls += 2*j

		qc.measure(q, res)

		backend = Aer.get_backend('aer_simulator')

		transpiled = transpile(qc, backend = backend)
		counts = backend.run(transpiled, shots = 1).result().get_counts()

		for val, count in counts.items():
			if val == '1':
				return no_calls


