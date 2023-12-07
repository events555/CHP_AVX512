import re
class Qudit:
    def __init__(self, dimension, a=None, b=None):
        self.dimension = dimension
        self.pauli_string = "Z" # assume computational basis
        if a is None or b is None:
            self.a, self.b = self.find_ab(dimension)
        else:
            self.a = a
            self.b = b
    def set(self, pauli_string):
        self.pauli_string = pauli_string
    def apply(self, gate, target=None):
        if gate == "SUM" and target is not None:
            target.pauli_string += re.findall('X', self.pauli_string)
            self.pauli_string += 'Z' * (self.dimension - 1) * len(re.findall('Z', target.pauli_string))
            return

        x_replacements = {'R': 'Z', 'P': 'XZ', 'S': 'X' * self.a}
        z_replacements = {'R': 'X' * (self.dimension - 1), 'P': 'Z', 'S': 'Z' * self.b}

        self.pauli_string = re.sub('X', lambda match: x_replacements[gate], self.pauli_string)
        self.pauli_string = re.sub('Z', lambda match: z_replacements[gate], self.pauli_string)

    def find_ab(self, d):
        if d == 2:
            return 1, 1
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

    def set(self, pauli_string):
        for i, qudit in enumerate(self.qudits):
            qudit.set(pauli_string[i])

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
        self.qudit_target_index = target


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
        if self.qudit_register.dimension==2:
            self.qudit_register = chp.simulate(self.qudit_register, self.gates) 
        else:
            for gate in self.gates:
                if gate.name == "SUM":
                    qudit = self.qudit_register[gate.qudit_index]
                    target = self.qudit_register[gate.qudit_target_index]
                    qudit.apply(gate.name, target)
                    target.pauli_string = self.simplify_string(target.pauli_string, target.dimension)
                else:
                    qudit = self.qudit_register[gate.qudit_index]
                    qudit.apply(gate.name)
                qudit.pauli_string = self.simplify_string(qudit.pauli_string, qudit.dimension)
    def __str__(self):
        gate_list = "\n".join([f"Gate: {gate.name}, Qudit Index: {gate.qudit_index}, Target Index: {gate.qudit_target_index}" for gate in self.gates])
        return f"Circuit:\n{gate_list}"
    def simplify_string(self, s, d):
        s = re.sub('X'*d, '' if len(s) > d else 'I', s)
        s = re.sub('Z'*d, '' if len(s) > d else 'I', s)
        if 'X' in s or 'Z' in s:
            s = s.replace('I', '')
        return s

qr = QuditRegister("Input Array", 2, 3)
qc = Circuit(qr)
qc.add_gate("R", 1)
qc.add_gate("P", 1)
qc.add_gate("P", 1)
qc.add_gate("SUM", 0, 1)
qc.simulate()
print(qr)