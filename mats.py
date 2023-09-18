from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product

pauli_matrices =                                                        \
        [                                                               \
                Matrix([[0, 1], [1, 0]]),                               \
                Matrix([[0,-I], [I, 0]]),                               \
                Matrix([[1, 0], [0, -1]])                               \
        ]

gell_mann_matrices =                                                    \
        [                                                               \
            Matrix([ [0, 1, 0], [1, 0, 0], [0, 0, 0] ]),                \
            Matrix([ [0,-I, 0], [I, 1, 0], [0, 0, 0] ]),                \
            Matrix([ [1, 0, 0], [0,-1, 0], [0, 0, 0] ]),                \
            Matrix([ [0, 0, 1], [0, 0, 0], [1, 0, 0] ]),                \
            Matrix([ [0, 0,-I], [0, 0, 0], [I, 0, 0] ]),                \
            Matrix([ [0, 0, 0], [0, 0, 1], [0, 1, 0] ]),                \
            Matrix([ [0, 0, 0], [0, 0,-I], [0,I , 0] ]),                \
            (1 / sqrt(3)) * Matrix([ [1+I*sqrt(3),   I*sqrt(3),   I*sqrt(3)],           \
                                      [-I*sqrt(3),   -2+I*sqrt(3), -2*I*sqrt(3)],       \
                                      [-I*sqrt(3), -2*I*sqrt(3),   -2+I*sqrt(3)] ]) 
        ]

def dft_gate(d):
    omega = exp((-2 * pi * I) / d)
    W = [[omega**(i * j) / sqrt(d) for j in range(d)] for i in range(d)]
    return Matrix(W)

def p_gate(d):
    omega = exp(pi*I/d)
    if d % 2 == 0:
        W = [[omega**((-i * (i + 2))) if i == j else 0 for j in range(d)] for i in range(d)]
    else:
        W = [[omega**((-i * (i + 1))) if i == j else 0 for j in range(d)] for i in range(d)]
    return Matrix(W)

def sqrt_Z(d):
    return sqrt(p_gate(d))

def sqrt_sqrt_Z(d):
    return sqrt(sqrt_Z(d))

#implement later
def SUM(d):
    if d == 2:
        return Matrix([[1, 0, 0 , 0], [0, 1, 0, 0],  [0, 0, 0 ,1], [0, 0, 1, 0]])
    return Matrix(d, d, lambda i, j: 1 if (i - j) % d == 0 else 0)

def chp(d, n):
    if n == 1:
        return [dft_gate(d), p_gate(d), eye(d)]
    elif n == 2:
        return [TensorProduct(*combination) for combination in product(chp(d,n-1), repeat=n)] + [SUM(d)]
    else:
        single_qudit_gates = chp(d, 1)
        gates = chp(d, n-1)
        return generate_matrices(single_qudit_gates, gates)

def generate_matrices(gates1, gates2):
    matrices = []
    for gate1 in gates1:
        for gate2 in gates2:
            matrix = TensorProduct(gate1, gate2)
            matrices.append(matrix)
    return matrices

if __name__ == "__main__":
    pprint(TensorProduct(SUM(2), eye(2)))
