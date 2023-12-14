import re

class QuditRegister:
    def __init__(self, name, dimension, num_qudits):
        self.name = name
        self.dimension = dimension
        self.num_qudits = num_qudits
        self.pauli_product = ['I' for _ in range(num_qudits)]

    def __getitem__(self, index):
        return self.pauli_product[index]

    def __str__(self):
            tensor_product = "\u2297 ".join(self.pauli_product)
            return f"{self.name}: {tensor_product}"

    def set(self, pauli, index=None):
        if index is None:
            pauli_split = pauli.split()
            for i in range(self.num_qudits):
                self.pauli_product[i] = pauli_split[i]
                simplified = self.simplify_string(self.pauli_product[i], self.dimension)
                if simplified is not None:
                    self.pauli_product = simplified
        else:
            self.pauli_product[index] = pauli
            simplified = self.simplify_string(self.pauli_product[index], self.dimension)
            if simplified is not None:
                self.pauli_product = simplified


        
    def simplify_string(self, s, d):
        s = re.sub('X'*d, '' if len(s) > d else 'I', s)
        s = re.sub('Z'*d, '' if len(s) > d else 'I', s)
        if 'X' in s or 'Z' in s:
            s = s.replace('I', '')

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
        if gate_name not in ["R", "P", "SUM"]:
            raise ValueError(f"Invalid gate: {gate_name}")
        elif gate_name == "SUM" and target_index not in range(self.qudit_register.num_qudits):
            raise ValueError(f"Invalid target index: {target_index}")
        elif qudit_index not in range(self.qudit_register.num_qudits):
            raise ValueError(f"Invalid qudit index: {qudit_index}")
        gate = Gate(gate_name, qudit_index, target_index)
        self.gates.append(gate)
    def __str__(self):
        gate_list = "\n".join([f"Gate: {gate.name}, Qudit Index: {gate.qudit_index}, Target Index: {gate.qudit_target_index}" for gate in self.gates])
        return f"Circuit:\n{gate_list}"
    