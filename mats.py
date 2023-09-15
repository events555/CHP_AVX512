import numpy as np
import math

imag_unit = complex(0+1j)

pauli_matrices =                                                        \
        [                                                               \
                np.array([[0+0j, 1+0j], [1+0j, 0+0j]]),                 \
                np.array([[0,0-1j], [0+1j, 0]]),                        \
                np.array([[1+0j, 0+0j], [0+0j, -1+0j]])                 \
        ]


gell_mann_matrices = \
        [
            np.array([ [0+0j, 1+0j, 0+0j], [1+0j, 0+0j, 0+0j], [0+0j, 0+0j, 0+0j] ]),           \
            np.array([ [0+0j, 0-1j, 0+0j], [1+0j, 0+1j, 0+0j], [0+0j, 0+0j, 0+0j] ]),           \
            np.array([ [1+0j, 0+0j, 0+0j], [0+0j, -1+0j, 0+0j], [0+0j, 0+0j, 0+0j] ]),          \
            np.array([ [0+0j, 0+0j, 1+0j], [0+0j, 0+0j, 0+0j], [1+0j, 0+0j, 0+0j] ]),           \
            np.array([ [0+0j, 0+0j, 0-1j], [0+0j, 0+0j, 0+0j], [0+1j, 0+0j, 0+0j] ]),           \
            np.array([ [0+0j, 0+0j, 0+0j], [0+0j, 0+0j, 1+0j], [0+0j, 1+0j, 0+0j] ]),           \
            np.array([ [0+0j, 0+0j, 0+0j], [0+0j, 0+0j, 0-1j], [0+0j, 0+1j, 0+0j] ]),           \
            (1 / math.sqrt(3)) * np.array([ [1+0j, 0+0j, 0+0j], [0+0j, 1+0j, 0+0j], [0+0j, 0+0j, -2+0j] ]),          \
        ]


def dft_gate(d):
    omega = np.exp((- 2 * math.pi * imag_unit) / d )
    W = np.eye(d, d, dtype='cfloat')
    it = np.nditer(W, flags = ['multi_index'], op_flags = ['readwrite'])
    for k in it:
        i = it.multi_index[0]
        j = it.multi_index[1]
        W[i, j] = np.power(omega, i * j) / math.sqrt(d)
        
    return W


def p_gate(d):
    omega = np.exp((- 2 * math.pi * imag_unit) / d )
    W = np.eye(d, d, dtype='cfloat')
    for k in range(d):
        W[k, k] = np.power(omega, ( k * (k - 1) ) // 2)

    return W

def sqrt_Z(d):
    return np.diag( [ np.sqrt( np.exp(( -2 * math.pi * imag_unit * k) / d) ) for k in range(d) ] )

def sqrt_sqrt_Z(d):
    return np.sqrt(sqrt_Z(d))

#implement later
def SUM(d):
   return np.eye(d, d, dtype='cfloat')

def chp(d):
    return [SUM(d), dft_gate(d), sqrt_Z(d)]

def chp_plus_labels(d):
    return [chp(d), ["C", "D", "Z"]]


qubit_comp_basis =                                                      \
        [                                                               \
            np.array([[0+0j], [1+0j]]),                                 \
            np.array([[1+0j], [0+0j]])                                  \
        ]

