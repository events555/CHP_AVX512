import mats
from sympy import *

def qudit_stab(qudits, gates, rounds):
    vectors = {q: q for q in qudits}
    for _ in range(rounds):
        for p in gates:
            for v in list(vectors.keys()):
                vec = p * v
                relPhase_vec = discard_global_phase_state(vec)
                if relPhase_vec not in vectors:
                    vectors[relPhase_vec] = relPhase_vec

    return vectors

# Function to discard global phase of a state
def discard_global_phase_state(mat):
    mat = mat.as_mutable()  # Convert to mutable matrix
    # Find the first non-zero element
    for i in range(mat.rows):
        if mat[i] != 0:  # Avoid division by zero
            global_phase = arg(mat[i])
            break
    # Discard the global phase from all elements
    for i in range(mat.rows):
        if mat[i] != 0:  # Avoid division by zero
            mat[i] = mat[i] * exp(-I * global_phase)
    return mat.as_immutable()  # Convert back to immutable matrix if needed

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

    while rounds_needed < 2:
        stab = qudit_stab(stab, gates, 1)
        rounds_needed += 1
    for keys in stab:
        print(keys)
    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
