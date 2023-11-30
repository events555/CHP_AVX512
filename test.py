from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product
from circuit import Qudit, QuditRegister, Gate, Circuit

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

import random

def statevec_test(d, num_qudits=2, num_gates=2):
    qr = QuditRegister("Test", d, num_qudits)
    qc = Circuit(qr)
    X, Z = generate_paulis(d)
    P, R, SUM, S = generate_clifford(d)

    # Initialize the statevector
    statevec = Matrix([1 if i == 0 else 0 for i in range(d**num_qudits)])

    # Add a series of random gates to the quantum circuit
    for _ in range(num_gates):
        gate_name = random.choice(["R", "P"] if num_qudits >= 2 else ["R", "P"])
        qudit_index = random.randrange(num_qudits)
        if gate_name == "SUM":
            qc.add_gate(gate_name, num_qudits-2, num_qudits-1)
        else:
            qc.add_gate(gate_name, qudit_index)

        # Generate the gate matrix
        if gate_name == "R":
            gate_matrix = R
        elif gate_name == "P":
            gate_matrix = P
        elif gate_name == "SUM":
            gate_matrix = SUM
        elif gate_name == "S":
            gate_matrix = S

        # Generate the list of matrices for the tensor product
        if gate_name == "SUM":
            matrices = [eye(d) for _ in range(num_qudits - 2)] + [gate_matrix, gate_matrix]
        else:
            matrices = [gate_matrix if i == qudit_index else eye(d) for i in range(num_qudits)]
        # Calculate the tensor product of the matrices
        tensor_product = matrices[0]
        for matrix in matrices[1:]:
            tensor_product = TensorProduct(tensor_product, matrix)
        # Apply the gate to the statevector
        statevec = tensor_product * statevec

     # Check if the final stabilizer stabilizes the statevector
    qc.simulate()
    stabilizers = []
    for i in range(qr.size):
        pauli_string = qr[i].pauli_string
        stabilizer = eye(d)
        for char in pauli_string:
            if char == 'X':
                stabilizer = stabilizer * X
            elif char == 'Z':
                stabilizer = stabilizer * Z
        stabilizers.append(stabilizer)

    final_stabilizer = stabilizers[0]
    for stabilizer in stabilizers[1:]:
        final_stabilizer = TensorProduct(final_stabilizer, stabilizer)

    if final_stabilizer * statevec != statevec:
        print(f"Not stabilized.")
    else:
        print(f"Stabilized.")
    pprint(final_stabilizer)
    pprint(statevec)
    print(qc)
    print(qr)


statevec_test(2, 3, 4)


