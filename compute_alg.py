import mats
from sympy import *
from sympy.physics.quantum import TensorProduct
from sympy.physics.quantum import tensor_product_simp

def qudit_stab(qudits, gates, rounds):
    # Convert qudits to numerical form and use tuple to make them hashable
    vectors = {q.evalf(): q for q in qudits}
    for _ in range(rounds):
        for p in gates:
            for v in list(vectors.keys()):
                vec = p * v
                relPhase_vec = discard_global_phase_state(vec)
                relPhase_vec_num = relPhase_vec.evalf()
                relPhase_vec_num_chopped = relPhase_vec_num.applyfunc(chop_small_numbers)
                if relPhase_vec_num_chopped not in vectors:
                    vectors[relPhase_vec_num_chopped] = vec
    return vectors

def chop_small_numbers(expr, epsilon=1e-10):
    re = expr.as_real_imag()[0]
    im = expr.as_real_imag()[1]

    if abs(re) < epsilon:
        re = 0
    if abs(im) < epsilon:
        im = 0

    return re + im*I

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


def n_qudit_comp_basis(d, n):
    single_qudit_basis = [ImmutableMatrix(d, 1, lambda i, j: 1 if i == k else 0) for k in range(d)]
    return [TensorProduct(*basis_vector) for basis_vector in product(single_qudit_basis, repeat=n)]

def two_qubit_comp_basis():
    # Define the single-qubit basis vectors
    single_qubit_basis = [ImmutableMatrix([[1], [0]]), ImmutableMatrix([[0], [1]])]
    # Generate all combinations of basis vectors
    return [TensorProduct(b1.as_mutable(), b2.as_mutable()).as_immutable() for b1 in single_qubit_basis for b2 in single_qubit_basis]

if __name__ == "__main__":
    
    #sanitize input later
    d = int(input("Dim: "))
    # n = int(input("Num: "))
    #gates = mats.pauli_matrices if d == 2 else mats.chp_3[:7]
    #initialization 
    gates = mats.chp(d, 2)
    #stab = n_qudit_comp_basis(d,n)
    stab = two_qubit_comp_basis()
    rounds_needed = 0

    while rounds_needed < 2:
        stab = qudit_stab(stab, gates, 1)
        rounds_needed += 1
    for key, value in stab.items():
        print(key)
    print(len(stab))
    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
