from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product


def generate_paulis(d):

    # Shift Matrix
    X = Matrix.zeros(d)
    for i in range(d):
        X[i, (i-1)%d] = 1

    # Clock Matrix
    Z = Matrix.eye(d)
    for i in range(d):
        Z[i, i] = exp(2*pi*I*i/d)

    return X, Z

def generate_clifford(d):
    # P Gate
    P = Matrix.eye(d)
    for j in range(d):
        P[j, j] = exp(I*pi*j*(j-1)/d)

    # DFT Matrix
    W = Matrix.zeros(d)
    for m in range(d):
        for n in range(d):
            W[m, n] = 1/sqrt(d) * exp(2*pi*I*m*n/d)

    # SUM Gate
    SUM = Matrix.zeros(d**2)
    for i, j in product(range(d), repeat=2):
        SUM[d*i + j, d*i + (i+j)%d] = 1
    SUM = SUM.reshape(d**2, d**2)

    # S Gate
    S = Matrix.zeros(d)
    for a in range(d):
        for b in range(d):
            if a*b % d == 1:
                S[a, a] = 1
                break

    return P, W, SUM, S
def generate_powers(X, Z, W, d):
    X_list = [(X**i).applyfunc(nsimplify) for i in range(1, d)]
    Z_list = [(Z**i).applyfunc(nsimplify) for i in range(1, d)]
    W_list = [(W**i).applyfunc(nsimplify) for i in range(1, d)]

    return X_list, Z_list, W_list

# Usage
d = 3  # Dimension
X, Z = generate_paulis(d)
P, W, SUM, S = generate_clifford(d)
X_list, Z_list, W_list = generate_powers(X, Z, W, d)
XZ_list = [(X_list[i-1] * Z_list[j-1]).applyfunc(nsimplify) for i in range(1, d) for j in range(1, d)]
def generate_W_mapping(X_list, Z_list, W):
    mapping = {}
    for gate_list, gate_name in [(X_list, 'X'), (Z_list, 'Z')]:
        for i, gate in enumerate(gate_list):
            conjugated_gate = W * gate * W.inv()
            conjugated_gate = conjugated_gate.applyfunc(nsimplify)
            for j, X_prime in enumerate(X_list):
                if conjugated_gate == X_prime:
                    mapping[f"{gate_name}^{i+1}"] = f"X^{j+1}"
                    break
            for j, Z_prime in enumerate(Z_list):
                if conjugated_gate == Z_prime:
                    mapping[f"{gate_name}^{i+1}"] = f"Z^{j+1}"
                    break
    return mapping
def generate_P_mapping(X_list, Z_list, XZ_list, P):
    mapping = {}
    for gate_list, gate_name in [(X_list, 'X'), (Z_list, 'Z')]:
        for i, gate in enumerate(gate_list):
            conjugated_gate = P * gate * P.inv()
            conjugated_gate = conjugated_gate.applyfunc(nsimplify)
            for j, X_prime in enumerate(X_list):
                if conjugated_gate == X_prime:
                    mapping[f"{gate_name}^{i+1}"] = f"X^{j+1}"
                    break
            for j, Z_prime in enumerate(Z_list):
                if conjugated_gate == Z_prime:
                    mapping[f"{gate_name}^{i+1}"] = f"Z^{j+1}"
                    break
            for j, (k, l) in enumerate(product(range(1, d), repeat=2)):
                if conjugated_gate == XZ_list[j]:
                    mapping[f"{gate_name}^{i+1}"] = f"X^{k}*Z^{l}"
                    break
    return mapping
def generate_S_mapping(X_list, Z_list, XZ_list, S):
    mapping = {}
    for gate_list, gate_name in [(X_list, 'X'), (Z_list, 'Z')]:
        for i, gate in enumerate(gate_list):
            conjugated_gate = S * gate * S.inv()
            conjugated_gate = conjugated_gate.applyfunc(nsimplify)
            for j, X_prime in enumerate(X_list):
                if conjugated_gate == X_prime:
                    mapping[f"{gate_name}^{i+1}"] = f"X^{j+1}"
                    break
            for j, Z_prime in enumerate(Z_list):
                if conjugated_gate == Z_prime:
                    mapping[f"{gate_name}^{i+1}"] = f"Z^{j+1}"
                    break
            for j, (k, l) in enumerate(product(range(1, d), repeat=2)):
                if conjugated_gate == XZ_list[j]:
                    mapping[f"{gate_name}^{i+1}"] = f"X^{k}*Z^{l}"
                    break
    return mapping
# Usage
W_mapping = generate_W_mapping(X_list, Z_list, W)
for key, value in W_mapping.items():
    print(f"{key} -> {value}")

print("\nP mapping:")
P_mapping = generate_P_mapping(X_list, Z_list, XZ_list, P)
for key, value in P_mapping.items():
    print(f"{key} -> {value}")