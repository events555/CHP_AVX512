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
                if relPhase_vec not in stabilzed:
                    stabilzed[relPhase_vec] = relPhase_vec
                    if len(stabilzed) % 10 == 0:
                        print(str(len(stabilzed)) + " states found.")
                        print(relPhase_vec)
                        print(vec)
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

    test3 = (1/4 - I/4)*exp(I*pi/4) #sqrt2/4
    test4 = sqrt(2)/4
    print(test3.equals(test4))
    print(simplify(test3-test4) == 0 )
    print(test3.evalf() == test4.evalf())
    print(simplify(arg((1 + I))))
    rounds_needed = 0
    while rounds_needed < 3:
        stab = qudit_stab(stab, gates, 1)
        rounds_needed += 1
    print(len(stab))
    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
