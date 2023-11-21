class Qudit:
    def __init__(self, dimension):
        self.dimension = dimension
        self.pauli_string = "X" # assume computational basis
    def apply(self, gate, target=None):
        if gate == "R":
            self.apply_R_gate()
        elif gate == "P":
            self.apply_P_gate()
        elif gate == "SUM" and target is not None:
            self.apply_SUM_gate(target)
    def apply_R_gate(self):
        if self.pauli_string == "X":
            self.pauli_string = "Z"
        elif self.pauli_string == "Z":
            self.pauli_string = "X^-1"
    def apply_P_gate(self):
        if self.pauli_string == "X":
            self.pauli_string = "XZ"
    def apply_SUM_gate(self, control_index, target_index):
        control_qudit = self.qudits[control_index]
        target_qudit = self.qudits[target_index]
        if control_qudit.pauli_string == "X" and target_qudit.pauli_string == "I":
            target_qudit.pauli_string = "X"
        elif control_qudit.pauli_string == "I" and target_qudit.pauli_string == "Z":
            control_qudit.pauli_string = "Z^-1"

class QuditRegister:
    def __init__(self, name, dimension, size):
        self.name = name
        self.dimension = dimension
        self.size = size
        self.qudits = [Qudit(dimension) for _ in range(size)]

    def __getitem__(self, index):
        return self.qudits[index]

    def __str__(self):
            tensor_product = "\u2297 ".join([qudit.pauli_string for qudit in self.qudits])
            return f"{self.name}: {tensor_product}"

class Gate:
    def __init__(self, name, qudit, target=None):
        self.name = name
        self.qudit_index = qudit
        self.qudit_target = target

class Circuit:
    def __init__(self, qudit_register):
        self.qudit_register = qudit_register
        self.gates = []

    def add_gate(self, gate_name, qudit_index, target_index=None):
        if gate_name not in ["R", "P", "SUM"]:
            raise ValueError(f"Invalid gate: {gate_name}")
        gate = Gate(gate_name, qudit_index, target_index)
        self.gates.append(gate)

    def simulate(self):
        for gate in self.gates:
            qudit = self.qudit_register[gate.qudit_index]
            if gate.name == "R" or gate.name == "P":
                qudit.apply(gate.name)
            elif gate.name == "SUM":
                target = self.qudit_register[gate.target_index]
                qudit.apply(gate.name, target)



qr = QuditRegister("Input Array", 3, 3)
qc = Circuit(qr)
qc.add_gate("R", 0)
qc.add_gate("P", 1)
print(qr)
qc.simulate()
print(qr)