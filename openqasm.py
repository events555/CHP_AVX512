import re

class Qudit:
    def __init__(self, dimension, a=None, b=None):
        self.dimension = dimension
        self.pauli_string = "X" # assume computational basis
        if a is None or b is None:
            self.a, self.b = self.find_ab(dimension)
        else:
            self.a = a
            self.b = b

    def apply(self, gate, target=None):
        new_pauli_string = ""
        for char in self.pauli_string:
            if char == "X":
                if gate == "R":
                    new_pauli_string += "Z"
                elif gate == "P":
                    new_pauli_string += "XZ"
                elif gate == "S":
                    new_pauli_string += "X"*self.a
            elif char == "Z":
                if gate == "R":
                    new_pauli_string += "X"*(self.dimension-1)
                elif gate == "S":
                    new_pauli_string += "Z"*self.b
        self.pauli_string = new_pauli_string
        if gate == "SUM" and target is not None:
            self.apply_SUM_gate(target)

    def apply_SUM_gate(self, control_index, target_index):
        control_qudit = self.qudits[control_index]
        target_qudit = self.qudits[target_index]
        if control_qudit.pauli_string == "X" and target_qudit.pauli_string == "I":
            target_qudit.pauli_string = "X"
        elif control_qudit.pauli_string == "I" and target_qudit.pauli_string == "Z":
            control_qudit.pauli_string = "Z"*(self.dimension-1)

    def find_ab(self, d):
        for a in range(1, d):
            for b in range(1, d):
                if (a * b) % d == 1 and not (a == 1 and b == 1):
                    return a, b
        return None, None

class QuditRegister:
    def __init__(self, name, dimension, size):
        self.name = name
        self.dimension = dimension
        self.size = size
        a, b = self.find_ab(dimension)
        self.qudits = [Qudit(dimension, a, b) for _ in range(size)]

    def __getitem__(self, index):
        return self.qudits[index]

    def __str__(self):
            tensor_product = "\u2297 ".join([qudit.pauli_string for qudit in self.qudits])
            return f"{self.name}: {tensor_product}"
    def find_ab(self, d):
        for a in range(1, d):
            for b in range(1, d):
                if (a * b) % d == 1 and not (a == 1 and b == 1):
                    return a, b
        return None, None

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
        if gate_name not in ["R", "P", "SUM", "S"]:
            raise ValueError(f"Invalid gate: {gate_name}")
        gate = Gate(gate_name, qudit_index, target_index)
        self.gates.append(gate)

    def simulate(self):
        for gate in self.gates:
            qudit = self.qudit_register[gate.qudit_index]
            if gate.name == "SUM":
                target = self.qudit_register[gate.target_index]
                qudit.apply(gate.name, target)
            else:
                qudit.apply(gate.name)
            qudit.pauli_string = self.replace_consecutive_chars(qudit.pauli_string, qudit.dimension)

    def replace_consecutive_chars(self, s, d):
        s = re.sub('X'*d, '' if len(s) > d else 'I', s)
        s = re.sub('Z'*d, '' if len(s) > d else 'I', s)
        return s

qr = QuditRegister("Input Array", 3, 3)
qc = Circuit(qr)
qc.add_gate("R", 0)
qc.add_gate("P", 1)
qc.add_gate("S", 2)
qc.add_gate("S", 0)
qc.add_gate("S", 0)
qc.add_gate("S", 2)
qc.add_gate("R", 0)
qc.add_gate("S", 0)
print(qr)
qc.simulate()
print(qr)