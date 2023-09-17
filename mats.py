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
def SUM(d, n):
    if d == 2:
        return Matrix([[1, 0, 0 , 0], [0, 1, 0, 0],  [0, 0, 0 ,1], [0, 0, 1, 0]])
    return Matrix(d, d, lambda i, j: 1 if (i - j) % d == 0 else 0)

def chp(d, n):
    gates = [dft_gate(d), p_gate(d), eye(d)]
    chp = n_qubit_gates(gates, n)
    if n > 1:
        chp.append(SUM(d, n))
    return chp


def n_qubit_gates(gates, n):
    if n == 1:
        return gates
    tensor_products = [TensorProduct(*combination) for combination in product(gates, repeat=n)]
    return tensor_products
