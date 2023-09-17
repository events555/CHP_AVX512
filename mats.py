from sympy import *
from sympy.physics.quantum import TensorProduct

pauli_matrices =                                                        \
        [                                                               \
                ImmutableMatrix([[0, 1], [1, 0]]),                               \
                ImmutableMatrix([[0,-I], [I, 0]]),                               \
                ImmutableMatrix([[1, 0], [0, -1]])                               \
        ]

gell_mann_matrices =                                                    \
        [                                                               \
            ImmutableMatrix([ [0, 1, 0], [1, 0, 0], [0, 0, 0] ]),                \
            ImmutableMatrix([ [0,-I, 0], [I, 1, 0], [0, 0, 0] ]),                \
            ImmutableMatrix([ [1, 0, 0], [0,-1, 0], [0, 0, 0] ]),                \
            ImmutableMatrix([ [0, 0, 1], [0, 0, 0], [1, 0, 0] ]),                \
            ImmutableMatrix([ [0, 0,-I], [0, 0, 0], [I, 0, 0] ]),                \
            ImmutableMatrix([ [0, 0, 0], [0, 0, 1], [0, 1, 0] ]),                \
            ImmutableMatrix([ [0, 0, 0], [0, 0,-I], [0,I , 0] ]),                \
            (1 / sqrt(3)) * ImmutableMatrix([ [1+I*sqrt(3),   I*sqrt(3),   I*sqrt(3)],           \
                                      [-I*sqrt(3),   -2+I*sqrt(3), -2*I*sqrt(3)],       \
                                      [-I*sqrt(3), -2*I*sqrt(3),   -2+I*sqrt(3)] ]) 
        ]

def dft_gate(d):
    omega = exp((-2 * pi * I) / d)
    W = [[omega**(i * j) / sqrt(d) for j in range(d)] for i in range(d)]
    return ImmutableMatrix(W)

def p_gate(d):
    omega = exp(pi*I/d)
    if d % 2 == 0:
        W = [[omega**((-i * (i + 2))) if i == j else 0 for j in range(d)] for i in range(d)]
    else:
        W = [[omega**((-i * (i + 1))) if i == j else 0 for j in range(d)] for i in range(d)]
    return ImmutableMatrix(W)

def sqrt_Z(d):
    return sqrt(p_gate(d))

def sqrt_sqrt_Z(d):
    return sqrt(sqrt_Z(d))

#implement later
def SUM(d):
    if d == 2:
        return ImmutableMatrix([[1, 0, 0 , 0], [0, 1, 0, 0],  [0, 0, 0 ,1], [0, 0, 1, 0]])
    return ImmutableMatrix(d, d, lambda i, j: 1 if (i - j) % d == 0 else 0)

def chp(d, n):
    chp = n_qubit_gates(n)
    chp.append(SUM(d))
    return chp


def n_qubit_gates(n):
    gates = [Matrix([[sqrt(2)/2, sqrt(2)/2], [sqrt(2)/2, -sqrt(2)/2]]), Matrix([[1,0],[0,I]]), Matrix([[1,0],[0,1]])]
    tensor_products = []
    for gate1 in gates:
        for gate2 in gates:
            tensor_products.append(TensorProduct(gate1, gate2).as_immutable())
    return tensor_products


qubit_comp_basis =                                                      \
        [                                                               \
            ImmutableMatrix([[0], [1]]),                                         \
            ImmutableMatrix([[1], [0]])                                          \
        ]