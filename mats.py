from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product

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
    omega = exp(pi*I/(2*d))
    if d % 2 == 0:
        W = [[omega**((-i * (i + 2))) if i == j else 0 for j in range(d)] for i in range(d)]
    else:
        W = [[omega**((-i * (i + 1))) if i == j else 0 for j in range(d)] for i in range(d)]
    return Matrix(W)

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
        return [TensorProduct(*combination) for combination in product(chp(d,1), repeat=2)] + [SUM(d)]
    else:
        single_qubit_gates = chp(d, 1)
        gates = chp(d, n-1)
        return [TensorProduct(*combination) for combination in product(single_qubit_gates, gates)]

def chpt(d, n):
    if n == 1:
        return [dft_gate(d), p_gate(d), eye(d), sqrt_Z(d)]
    elif n == 2:
        return [TensorProduct(*combination) for combination in product(chpt(d,1), repeat=2)] + [SUM(d)]
    else:
        single_qubit_gates = chpt(d, 1)
        gates = chpt(d, n-1)
        return [TensorProduct(*combination) for combination in product(single_qubit_gates, gates)]

if __name__ == "__main__":
    pprint(chpt(2,2))
