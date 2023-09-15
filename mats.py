from sympy import *

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
    omega = exp((-2 * pi * I) / d)
    W = [[omega**((k * (k - 1)) // 2) if i == j else 0 for j in range(d)] for i in range(d)]
    return ImmutableMatrix(W)

def sqrt_Z(d):
    return ImmutableMatrix.diag(*[sqrt(exp((-2 * pi * I * k) / d)) for k in range(d)])

# might need to convert to mutable matrix and then back to immutable matrix
def sqrt_sqrt_Z(d):
    return sqrt(sqrt_Z(d))

#implement later
def SUM(d):
    return ImmutableMatrix(d, d, lambda i, j: 1 if (i - j) % d == 0 else 0)

def chp(d):
    return [SUM(d), dft_gate(d), sqrt_Z(d)]

def chp_plus_labels(d):
    return [chp(d), ["C", "D", "Z"]]


qubit_comp_basis =                                                      \
        [                                                               \
            ImmutableMatrix([[0], [1]]),                                         \
            ImmutableMatrix([[1], [0]])                                          \
        ]

