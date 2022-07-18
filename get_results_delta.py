from abstract_quantum_simulation import *

H = 40
diff = math.sqrt(1/H)
means = [0.5 - diff/2, 0.5 + diff/2]
vals = []
for delta100 in range(3, 31, 3):
    delta = delta100/100
    print("#####\nEvaluating for delta =", delta)
    vals.append(BestArm(means, delta))

    print("vals so far:", vals)

# [42581415603, 38991732519, 38102359083, 35786837118, 35503310730, 34048756866, 32950261698, 32919858150, 32259092025, 31916263149]