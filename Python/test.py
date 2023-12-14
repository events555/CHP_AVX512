from sympy import *
from sympy.physics.quantum import TensorProduct
from itertools import product
from circuit import QuditRegister, Circuit
from tableau import Tableau, Program
import time
import random

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
    if d == 2:
        P[d-1, d-1] = exp(I*pi/2)
    else:
        P[d-1, d-1] = exp(I*2*pi/d)

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

    return P, R, SUM

def discard_global_phase_state(mat):
    mat = mat.as_mutable()  # Convert to mutable matrix
    global_phase = None
    for i in range(mat.rows):
        if N(mat[i], 15, chop=True) != 0:  # Avoid division by zero, values close to 0 were getting through != and arg was returning NaN for 0
            global_phase = arg(mat[i])
            break
    if global_phase is not None:
        # Discard the global phase from all elements
        for i in range(mat.rows):
            mat[i] = mat[i] * exp(-I * global_phase)
    mat_num = N(mat, 15, chop=True)
    mat_simplify = mat_num.applyfunc(nsimplify)
    return mat_simplify.as_immutable()  # Convert back to immutable matrix

def statevec_test(d, trials=1, num_qudits=2, num_gates=2, seed=int(time.time())):
    for trial in range(trials):
        print("Trial %d" % trial)
        # Initialize the quantum program
        qr = QuditRegister("Trial %d" % trial, d, num_qudits)
        qc = Circuit(qr)
        table = Tableau(d, num_qudits)
        # Generate the Pauli and Clifford gates
        random.seed(seed+trial)
        X, Z = generate_paulis(d)
        P, R, SUM = generate_clifford(d)
        # Initialize the statevector
        statevec = Matrix([1 if i == 0 else 0 for i in range(d**num_qudits)])
        # Add a series of random gates to the quantum circuit
        for _ in range(num_gates):
            gate_name = random.choice(["R", "P"] if num_qudits >= 2 else ["R", "P"])
            qudit_index = random.randrange(num_qudits-1)
            if gate_name == "SUM":
                target_index = random.choice([i for i in range(num_qudits) if i != qudit_index])
                qc.add_gate(gate_name, qudit_index, target_index)
            else:
                qc.add_gate(gate_name, qudit_index)
            
            gate_map = {'R': R, 'P': P, 'SUM': SUM}

            # Generate the list of matrices for the tensor product
            if gate_name == "SUM":
                matrices = [eye(d) for _ in range(num_qudits - 2)] + [gate_map[gate_name], gate_map[gate_name]]
            else:
                matrices = [gate_map[gate_name] if i == qudit_index else eye(d) for i in range(num_qudits)]
            # Calculate the tensor product of the matrices
            tensor_product = matrices[0]
            for matrix in matrices[1:]:
                tensor_product = TensorProduct(tensor_product, matrix)
            # Apply the gate to the statevector
            statevec = (tensor_product * statevec).applyfunc(nsimplify)

        # Check if the final stabilizer stabilizes the statevector
        prog = Program(table, qc)
        prog.simulate()
        print(prog.stabilizer_tableau)
        print(qc)
        pauli_map = {'I': eye(d), 'X': X, 'Z': Z}
        pauli_matrices = []
        for j in range(num_qudits, 2*num_qudits):
            stabilizer = prog.get_stabilizer(j)
            for matrix in stabilizer.pauli_product:
                pauli_stabilizer = eye(d)
                for char in matrix:
                    pauli_stabilizer = pauli_stabilizer*pauli_map[char]
                pauli_matrices.append(pauli_stabilizer)
            pauli_stabilizer = pauli_matrices[0]
            for i in range(1, num_qudits):
                pauli_stabilizer = TensorProduct(pauli_stabilizer, pauli_matrices[i]).applyfunc(nsimplify)
            print("Stabilizer %d" % (j - num_qudits))
            pprint(statevec)
            pprint(pauli_stabilizer)
            stabilized = discard_global_phase_state(pauli_stabilizer * statevec)
            statevec = discard_global_phase_state(statevec)
            if not all(simplify(i) == 0 for i in (Matrix(stabilized) - Matrix(statevec))):
                raise Exception(f"Trial {trial + 1} was not stabilized. Seed: {seed}")



#statevec_test(2, 1, 4, 8, 1702577339)


