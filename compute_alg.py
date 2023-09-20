import mats
from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product

def qudit_stab(stabilized, gates, rounds):
    for v in list(stabilized.values()): # given every currently seen starting state, test all gates of depth round
        for _ in range(rounds):
            for p in gates:
                vec = p * v
                relPhase_vec = discard_global_phase_state(vec)
                if relPhase_vec not in stabilized:
                    stabilized[relPhase_vec] = relPhase_vec
                    if len(stabilized) % 10 == 0 and len(stabilized) >= 360:
                        print(relPhase_vec)
    return stabilized

def discard_global_phase_state(mat):
    mat = mat.as_mutable()  # Convert to mutable matrix
    global_phase = None
    for i in range(mat.rows):
        if N(mat[i], 15, chop=True) != 0:  # Avoid division by zero, values close to 0 were getting through != and arg was returning NaN for 0
            global_phase = arg(mat[i])
            break
    if global_phase is not None:
        # Discard the global phase from all elements
        for i in range(mat.rows):
            mat[i] = mat[i] * exp(-I * global_phase)
    mat_num = N(mat, 15, chop=True)
    mat_simplify = mat_num.applyfunc(nsimplify)
    return mat_simplify.as_immutable()  # Convert back to immutable matrix


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
    max_rounds = int(input("Max rounds: "))
    t_gate = input("Include T-gate? Y/N: ")
    #initialization 
    if t_gate.lower() == 'y':
        gates = mats.chpt(d, n)
    else:
        gates = mats.chp(d, n)
        
    qudits = n_qudit_comp_basis(n)
    for i in range(len(gates)):
        gates[i] = gates[i].as_immutable()
        #print(gates[i])
    for i in range(len(qudits)):
        qudits[i] = qudits[i].as_immutable()
        #print(stab[i])
    stabilized = {q: q for q in qudits}
    rounds_needed = 1
    prev_seen = len(qudits)
    while rounds_needed < max_rounds:
        stabilized = qudit_stab(stabilized, gates, rounds_needed)
        now_seen = len(stabilized)
        if prev_seen == now_seen:
            break
        else:
            prev_seen = now_seen
            rounds_needed += 1
            print("Round " + str(rounds_needed) + " complete. " + str(now_seen) + " states seen.")
    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
    # Ask the user if they want to print the seen states
    print_states = input("Would you like to print the discovered states? Y/N: ")
    
    if print_states.lower() == 'y':
        for key in stabilized:
            print(key.T)
