import math
import random
# from quantum_algo import *
from complex_subroutines import *


def algo_A_probs(means, l2, l1, alpha):
	"""
	Input:
	alpha: highest error allowed in the correctness of the estimate is 2*alpha/no_arms
	Output:
	no of calls to oracle AND
	list of (1-p1,p1) for each arm
	where p1 is probability of algorithm identifying the mean of the arm to be more than l1
	"""
	# for each arm return prob of getting 0 and 1 respectively
	no_arms = len(means)

	delta = l1 - l2
	m = math.ceil(math.log(1/delta, 2)) + 2

	a = alpha/(2*m*(no_arms**1.5))

	probs = []

	total_calls_to_oracle = 0

	for i, p in enumerate(means):

		# print("\n=====\nmean=", p)

		prob_of_1 = 0
		prob_of_0_sofar = 1

		for j in range(1, m + 1):
			eps = 2**-j
			M = 2*math.pi/(math.sqrt(1 + 2*eps) - 1)
			no_prec_qubits = math.ceil(math.log(M, 2))
			M = 2**no_prec_qubits

			if no_prec_qubits > 9:
				print("simulating total qubits=", no_prec_qubits + 1)

			no_calls, ae_sv = amplitude_est(
				ith_arm_oracle(i, means),
				1, no_prec_qubits, sim = 'statevector_simulator')


			# print("statevector:", list(enumerate(ae_sv)))

			prob_of_1_here = 0

			for val, ampl in enumerate(ae_sv):
				# divide val by 2 because of state qubit
				val = val//2
				estimate = math.sin(math.pi * val/(M))**2
				# estimate = math.sin(math.pi * val/(2*M))**2 # divide by 2M instead of M because of state qubit
				prob = abs(ampl)**2
				# print("estimate:", estimate, " with prob: ", prob)
				if estimate >= l1 - 3*eps/2:
					continue
				prob_of_1_here += prob

			# print("prob_of_1_here:", prob_of_1_here)

			ae_confidence = 8/(math.pi**2)
			## no of times to do ae to get good confidence on gae
			# no_ae = math.ceil(
			# 	math.log(2*a*a)/math.log(2*(1 - ae_confidence)))
			no_ae = math.ceil(math.log(2/a)/2*((ae_confidence - 1/2)**2))
			# print("no_ae: ", no_ae)

			no_calls *= no_ae
			total_calls_to_oracle += no_calls

			# calc prob of getting 1 if majority vote were taken
			majority_prob_1_here = 0
			for no_votes in range(no_ae//2 + 1, no_ae + 1):
				majority_prob_1_here += math.comb(no_ae, no_votes)*(
					prob_of_1_here**no_votes)*(
					(1 - prob_of_1_here)**(no_ae - no_votes))

			# print("majority_prob_1_here=", majority_prob_1_here)

			# print("eps=", eps, "prob_of_1_here:", prob_of_1_here)

			prob_of_1 += prob_of_0_sofar*majority_prob_1_here
			prob_of_0_sofar *= (1 - majority_prob_1_here)

		# print("overall prob of 1:", prob_of_1)


		# F is inverted
		probs.append((prob_of_1, 1 - prob_of_1))





		# if p < l1 - 2*2**(-m):
		# 	probs.append((1, 0))

		# else:
		# 	probs.append((0, 1))


	# print("==========")
	# print("(l2, l1)= ", l2, l1)
	# print("num calls by A=", total_calls_to_oracle)
	# print("l1 - eps = ", l1 - 2*2**(-m))
	# print(probs)
	# print("==========\n")
	return total_calls_to_oracle, probs

def estimate_probs(calls_probs, eps, delta):
	"""
	Input:
	calls_probs = 
	(no of calls to oracle by algo A, list of (1-p1, p1) for each arm. p1 is probability of success for some experiment on the arm)
	Output:
	total no of calls to oracle in estimating
	let est = sum_i(p1_i)/no_arms
	the algorithm returns estimate of est with atmost eps multiplicative error with probability >=1-delta
	"""
	calls_by_A, probs = calls_probs
	exact_prob = 0
	for i in probs:
		exact_prob += i[1]

	exact_prob /=len(probs)

	calls_to_A, estimated_prob = mult_bound_estimate(exact_prob, eps, -math.log(delta, 2))

	# print("estimate: no_calls_to_A=", calls_to_A)

	# print("exact:", exact_prob, "estimate:", estimated_prob)

	total_calls = calls_by_A*calls_to_A

	# return exact_prob
	return total_calls, estimated_prob

def amplify_probs(calls_probs):
	"""
	Input:
	calls_probs = 
	(no of calls to oracle by algo A, list of (1-p1, p1) for each arm. p1 is probability of success for some experiment on the arm)
	Output:
	let amp be a distribution over all arms where probability of each arm is proportional to its p1
	Only for VT version not implemented yet:
	# The procedure returns a sample from an estimate of amp
	# The estimated distribution is very accurate with probability >=1-delta
	"""
	calls_by_A, probs = calls_probs
	prob_of_1 = 0
	arm_selection_probs = []
	for prob in probs:
		arm_selection_probs.append(prob[1])
		prob_of_1 += prob[1]


	# valid_arms = []
	# for i in range(len(probs)):
	# 	# if probs is 1
	# 	if probs[i] == (0, 1):
	# 		valid_arms.append(i)

	print("Arm probs:", arm_selection_probs)
	calls_to_A = search(prob_of_1)
	total_calls = calls_by_A*calls_to_A

	return total_calls, random.choices(list(range(len(probs))), arm_selection_probs)[0]

	# print("Valid arms: ", valid_arms)
	# return random.choice(valid_arms)



def Shrink(means, k, I, delta):
	"""
	Input:
	means: means of all arms
	k: 1 or 2, depending on whether shrink is being called around highest or second highest arm
	I: interval to shrink
	Output:
	no of calls to oracle AND
	a subset interval of I that is 3/5 in length of I such that p_k lies in it with probability >=1-delta
	"""
	(a, b) = I
	n = len(means)
	eps = (b - a)/5
	delta = delta/2
	means = [1] + means

	calls1, r1 = estimate_probs(algo_A_probs(means, a + eps, a + 3*eps, 0.01*delta), 0.1, delta)
	calls2, r2 = estimate_probs(algo_A_probs(means, a + 2*eps, a + 4*eps, 0.01*delta), 0.1, delta)
	# print("Shrink: r1:", r1, " r2: ", r2)

	total_calls = calls1 + calls2

	# print("num_calls_in_shrink=", "{:e}".format(total_calls))

	B1 = 0
	B2 = 0

	if r1 > (k + 0.5)/(n + 1):
		B1 = 1
	if r2 > (k + 0.5)/(n + 1):
		B2 = 1

	if (B1, B2) == (0, 0):
		return total_calls, (a, a + 3*eps)
	elif (B1, B2) == (0, 1):
		return total_calls, (a + eps, a + 4*eps)
	elif (B1, B2) == (1, 0):
		return total_calls, (a + eps, a + 4*eps)
	else:
		return total_calls, (a + 2*eps, a + 5*eps)


def locate(means, delta):
	"""
	Output:
	no of calls to oracle AND
	two intervals i1, i2
	i1: highest mean lies in this interval
	i2: second highest mean in this interval
	min i1 -  max i2 >= 2|i1|
	"""
	i1 = (0, 1)
	i2 = (0, 1)
	delta = delta/8

	total_calls = 0

	num_calls_to_shrink = 0

	while i1[0] - i2[1] < 2*(i1[1] - i1[0]):
		calls1, i1 = Shrink(means, 1, i1, delta)
		calls2, i2 = Shrink(means, 2, i2, delta)

		num_calls_to_shrink += 2

		print("i1, i2: ", i1, i2)

		total_calls += (calls1 + calls2)

		delta = delta/2

	print("num_calls_to_shrink=", "{:e}".format(num_calls_to_shrink))

	return total_calls, i1, i2

def BestArm(means, delta):
	"""
	Input:
	means: list of means of arms
	Output:
	total number of calls to oracle
	index of arm with highest mean with probability >= 1-delta is printed
	"""
	delta = delta/2
	calls_for_locate, i1, i2 = locate(means, delta)

	print("total calls made for locate=", "{:e}".format(calls_for_locate))

	l1 = i1[0]
	l2 = i2[1]

	calls_for_search, index = amplify_probs(algo_A_probs(means, l2, l1, 0.01*delta))
	print("calls for search=", calls_for_search)
	calls = calls_for_search + calls_for_locate

	print("idx found=", index)

	return calls


