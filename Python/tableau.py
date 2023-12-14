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
        else:
            raise ValueError(f"Two qudit gates are not supported yet.")
class Program:
    def __init__(self, tableau, circuit):
        self.stabilizer_tableau = tableau
        self.circuit = circuit
    def simulate(self):
        for stabilizer in self.stabilizer_tableau.tableau:
            for gate in self.circuit.gates: 
                self.stabilizer_tableau.conjugate(gate.qudit_index, gate.qudit_target_index, gate.name, stabilizer)

table = Tableau(2, 2)
qudit_register = QuditRegister("Qudit Register", 2, 2)
circuit = Circuit(qudit_register)
circuit.add_gate("R", 0)
circuit.add_gate("P", 1)
#circuit.add_gate("SUM", 2, 3)
print(circuit)

program = Program(table, circuit)
program.simulate()
print(table)