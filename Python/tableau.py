from circuit import QuditRegister, Circuit
import re

class Tableau:
    def __init__(self, dimension, num_qudits):
        self.dimension = dimension
        self.num_qudits = num_qudits
        self.tableau = [QuditRegister("Destabilizer %d" % i, dimension, num_qudits) for i in range(num_qudits)]
        self.tableau.extend([QuditRegister("Stabilizer %d" % i, dimension, num_qudits) for i in range(num_qudits)])
        self.phase = [0 for _ in range(2*num_qudits)]
        for i in range(num_qudits):
            s = ['I'] * num_qudits  # Start with a list of 'I's
            s[i] = 'X'  # Replace the i-th element with 'X'
            self.tableau[i].set(' '.join(s))
            s[i] = 'Z'
            self.tableau[i + num_qudits].set(' '.join(s))

    def __str__(self):
        return "\n".join([f"{str(qudit_register)} Phase: {self.phase[i]}" for i, qudit_register in enumerate(self.tableau)])
    def __getitem__(self, index):
        return self.tableau[index]
    def conjugate(self, qudit_index, qudit_target_index, gate, stabilizer):
        # Conjugate's According to Gottesman's Higher Dimensional Mappings
        # https://arxiv.org/abs/quant-ph/9802007
        
        if qudit_target_index == None:
            if gate == "R":
                replacements = {'X': 'Z', 'Z': 'X' * (self.dimension - 1)}
                stabilizer.set(''.join([replacements.get(c, c) for c in stabilizer.pauli_product[qudit_index]]), qudit_index)
            elif gate == "P":
                replacements = {'X': 'XZ', 'Z': 'Z'}
                stabilizer.set(''.join([replacements.get(c, c) for c in stabilizer.pauli_product[qudit_index]]), qudit_index)
        elif qudit_target_index != None:
            if gate == "SUM":
                pass
                # control_string = list(stabilizer.pauli_product[qudit_index])
                # target_string = list(stabilizer.pauli_product[qudit_target_index])
                # for i, (c_control, c_target) in enumerate(zip(control_string, target_string)):
                #     if c_control == 'X':
                #         target_string[i] += 'X'
                #     if c_target == 'Z':
                #         control_string[i] += 'Z' * (self.dimension - 1)
                # stabilizer.set(''.join(control_string), qudit_index)
                # stabilizer.set(''.join(target_string), qudit_target_index)
        else:
            raise ValueError(f"No targets specified.")
    def check_commute(self, stab1, stab2):
        # Check if two stabilizers commute
        pauli1 = stab1.pauli_product
        pauli2 = stab2.pauli_product
        num_yx_or_zx_pairs = sum((p1 == 'XZ' or (p1 == 'X' and p2 == 'Z') or (p1 == 'Z' and p2 == 'X')) for p1, p2 in zip(pauli1, pauli2))
        return num_yx_or_zx_pairs % 2 == 0


class Program:
    def __init__(self, tableau, circuit):
        self.stabilizer_tableau = tableau
        self.num_qudits = tableau.num_qudits
        self.circuit = circuit
    def simulate(self):
        for gate in self.circuit.gates: 
            if gate.name == "MEASURE":
                for stabilizer in self.stabilizer_tableau.tableau[self.num_qudits:]:  # Only iterate over stabilizers
                    if stabilizer.pauli_product[gate.qudit_index] == 'Z':
                        stabilizer.set('I', gate.qudit_index)
            else:
                for stabilizer in self.stabilizer_tableau.tableau:
                    self.stabilizer_tableau.conjugate(gate.qudit_index, gate.qudit_target_index, gate.name, stabilizer)
            print(self.stabilizer_tableau)
            print("\n")

table = Tableau(2, 2)
qudit_register = QuditRegister("Qudit Register", 2, 2)
circuit = Circuit(qudit_register)
circuit.add_gate("R", 0)
circuit.add_gate("P", 0)
circuit.add_gate("SUM", 0, 1)
print(circuit) 

program = Program(table, circuit)
program.simulate()

#Revelation... page 3 of Gottesman paper says that XZ has order 2d for d=2 so it needs a factor of i, but d=odd does not