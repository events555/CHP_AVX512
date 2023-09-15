import mats
from sympy import *

def qudit_stab(qudits, gates, rounds):
    vectors = {q: q for q in qudits}
    for _ in range(rounds):
        for p in gates:
            for v in list(vectors.keys()):
                vec = p * v
                if vec not in vectors:
                    vectors[vec] = vec

    return vectors

def normalize_global_phase(matrix):
    max_element = max(matrix, key=abs)
    if max_element != 0:
        return simplify(matrix / max_element)
    return matrix

def qudit_comp_basis(d):
    return [ImmutableMatrix(d, 1, lambda i, j: 1 if i == k else 0) for k in range(d)]


if __name__ == "__main__":
    
    #sanitize input later
    d = int(input("Dim: "))
    
    #gates = mats.pauli_matrices if d == 2 else mats.chp_3[:7]
    
    #initialization 
    gates, symb = mats.chp_plus_labels(d)
    gate_seq = ["|" + str(k) + ">" for k in range(d)]
    prev_gates = 0
    stab = qudit_comp_basis(d)
    rounds_needed = 0

    while rounds_needed < 1:
        stab = qudit_stab(stab, gates, 1)
        rounds_needed += 1
    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
