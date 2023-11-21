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
        P[j, j] = exp(I*2*pi/d*j*(j-1)/2)

    # DFT Matrix
    R = Matrix.zeros(d)
    for m in range(d):
        for n in range(d):
            R[m, n] = 1/sqrt(d) * exp(2*pi*I*m*n/d)

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

    return P, R, SUM, S

def generate_powers(X, Z, d):
    X_list = [(X**i).applyfunc(nsimplify) for i in range(1, d)]
    Z_list = [(Z**i).applyfunc(nsimplify) for i in range(1, d)]

    return X_list, Z_list

def generate_tensor_products(X, Z, I):
    # Generate tensor products
    XI = TensorProduct(X, I)
    IX = TensorProduct(I, X)
    ZI = TensorProduct(Z, I)
    IZ = TensorProduct(I, Z)
    XX = TensorProduct(X, X)
    ZZ = TensorProduct(Z, Z)
    XZ = TensorProduct(X, Z)
    ZX = TensorProduct(Z, X)
    Z_InvZ = TensorProduct(Z.inv(), Z)
    return XI, IX, ZI, IZ, XX, ZZ, XZ, ZX, Z_InvZ

def generate_R_mapping(X_list, Z_list, R):
    mapping = {}
    R_inv = R.inv()
    for gate_list, gate_name in [(X_list, 'X'), (Z_list, 'Z')]:
        for i, gate in enumerate(gate_list):
            conjugated_gate = R * gate * R_inv
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
    P_inv = P.inv()
    for gate_list, gate_name in [(X_list, 'X'), (Z_list, 'Z')]:
        for i, gate in enumerate(gate_list):
            conjugated_gate = P * gate * P_inv
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
    S_inv = S.inv()
    for gate_list, gate_name in [(X_list, 'X'), (Z_list, 'Z')]:
        for i, gate in enumerate(gate_list):
            conjugated_gate = S * gate * S_inv
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

def generate_SUM_mapping(tensor_products, SUM):
    mapping = {}
    tensor_product_names = ['X ⊗ I', 'I ⊗ X', 'Z ⊗ I', 'I ⊗ Z', 'X ⊗ X', 'Z ⊗ Z', 'X ⊗ Z', 'Z ⊗ X', 'Z^-1 ⊗ Z']
    SUM_inv = SUM.inv()
    for i, tensor_product in enumerate(tensor_products):
        conjugated_gate = SUM * tensor_product * SUM_inv
        conjugated_gate = conjugated_gate.applyfunc(nsimplify)
        for j, tensor_product_prime in enumerate(tensor_products):
            if conjugated_gate.equals(tensor_product_prime):
                mapping[tensor_product_names[i]] = tensor_product_names[j]
                break
    return mapping

# Usage
d = 5  # Dimension
X, Z = generate_paulis(d)
P, R, SUM, S = generate_clifford(d)
X_list, Z_list = generate_powers(X, Z, d)
XZ_list = [(X_list[i-1] * Z_list[j-1]).applyfunc(nsimplify) for i in range(1, d) for j in range(1, d)]
tensor_products = generate_tensor_products(X, Z, Matrix.eye(d))

print("\nR mapping:")
R_mapping = generate_R_mapping(X_list, Z_list, R)
for key, value in R_mapping.items():
    print(f"{key} -> {value}")

print("\nP mapping:")
P_mapping = generate_P_mapping(X_list, Z_list, XZ_list, P)
for key, value in P_mapping.items():
    print(f"{key} -> {value}")


print("\nSUM mapping:")
SUM_mapping = generate_SUM_mapping(tensor_products, SUM)
for key, value in SUM_mapping.items():
    print(f"{key} -> {value}")

