from circuit import Qudit, QuditRegister, Gate, Circuit

class Tableau:
    def __init__(self, dimension, num_qudits):
        self.dimension = dimension
        self.num_qudits = num_qudits
        self.tableau = [QuditRegister("Stabilizer %d" % i, dimension, num_qudits) for i in range(num_qudits)]
        self.tableau.extend([QuditRegister("Destabilizer %d" % i, dimension, num_qudits) for i in range(num_qudits)])
        self.phase = [0 for _ in range(2*num_qudits)]
        for i in range(num_qudits):
            s = ['I'] * num_qudits  # Start with a list of 'I's
            s[i] = 'X'  # Replace the i-th element with 'X'
            self.tableau[i].set(''.join(s))
            s[i] = 'Z'
            self.tableau[i + num_qudits].set(''.join(s))

    def __str__(self):
        return "\n".join([f"{str(qudit_register)} Phase: {self.phase[i]}" for i, qudit_register in enumerate(self.tableau)])

    
    def __getitem__(self, index):
        return self.tableau[index]

table = Tableau(2, 4)
print(table)