import mats
from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product

def qudit_stab(qudits, gates, rounds):
    stabilzed = {}
    for _ in range(rounds):
        for p in gates:
            for v in qudits:
                vec = p * v
                relPhase_vec = discard_global_phase_state(vec)
                relPhase_vec_num = N(relPhase_vec, 15, chop=True)
                relPhase_vec_simplify = relPhase_vec_num.applyfunc(nsimplify)
                if relPhase_vec_simplify not in stabilzed:
                    stabilzed[relPhase_vec_simplify] = relPhase_vec_simplify
                    if len(stabilzed)% 80 == 0:
                        print(relPhase_vec_simplify)
    return stabilzed

def discard_global_phase_state(mat):
    mat = mat.as_mutable()  # Convert to mutable matrix
    global_phase = None
    for i in range(mat.rows):
        if mat[i] != 0:  # Avoid division by zero
            global_phase = arg(mat[i])
            break
    if global_phase is not None:
        # Discard the global phase from all elements
        for i in range(mat.rows):
            mat[i] = mat[i] * exp(-I * global_phase)
            mat[i].expand(Complex=True)
    return mat.as_immutable()  # Convert back to immutable matrix if needed


def n_qudit_comp_basis(n):
    # Define the single-qubit basis vectors
    single_qubit_basis = [Matrix(d, 1, lambda i, j: 1 if i == k else 0) for k in range(d)]
    if n == 1:
        return single_qubit_basis
    # Generate all combinations of basis vectors
    return [TensorProduct(*combination) for combination in product(single_qubit_basis, repeat=n)]


if __name__ == "__main__":
    
    #sanitize input later
    d = int(input("Dim: "))
    n = int(input("Num: "))
    #initialization 
    gates = mats.chp(d, n)
    stab = n_qudit_comp_basis(n)
    for i in range(len(gates)):
        gates[i] = gates[i].as_immutable()
        #print(gates[i])

    for i in range(len(stab)):
        stab[i] = stab[i].as_immutable()
        #print(stab[i])
    rounds_needed = 0
    while rounds_needed < 5:
        stab = qudit_stab(stab, gates, 1)
        rounds_needed += 1
        print("Round " + str(rounds_needed) + ": ")
    print(len(stab))
    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
