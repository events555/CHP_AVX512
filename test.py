from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product


def generate_matrices(d):
    # Clock Matrix
    U = Matrix.eye(d)
    for i in range(d):
        U[i, i] = exp(2*pi*I*i/d)

    # Shift Matrix
    V = Matrix.zeros(d)
    for i in range(d):
        V[i, (i+1)%d] = 1

    # DFT Matrix
    W = Matrix.zeros(d)
    for m in range(d):
        for n in range(d):
            W[m, n] = 1/sqrt(d) * exp(2*pi*I*m*n/d)

    return U, V, W

def generate_powers(U, V, W, d):
    U_list = [(U**i).applyfunc(nsimplify) for i in range(1, d)]
    V_list = [(V**i).applyfunc(nsimplify) for i in range(1, d)]
    W_list = [(W**i).applyfunc(nsimplify) for i in range(1, d)]

    return U_list, V_list, W_list

# Usage
d = 4  # Dimension
U, V, W = generate_matrices(d)
U_list, V_list, W_list = generate_powers(U, V, W, d)

# Print the matrices
for i, U in enumerate(U_list, start=1):
    print(f"U^{i}:")
    pprint(U)
    print()

for i, V in enumerate(V_list, start=1):
    print(f"V^{i}:")
    pprint(V)
    print()

for i, W in enumerate(W_list, start=1):
    print(f"W^{i}:")
    pprint(W)
    print()