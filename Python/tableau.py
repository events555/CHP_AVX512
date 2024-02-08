from paulistring import PauliString
from tableau_simulator import apply_gate as apply_gate


class Tableau:
    def __init__(self, num_qudits, dimension=2):
        self.num_qudits = num_qudits
        self.dimension = dimension
        self.xlogical = [
            PauliString(num_qudits, dimension=dimension) for _ in range(num_qudits)
        ]
        self.zlogical = [
            PauliString(num_qudits, dimension=dimension) for _ in range(num_qudits)
        ]

    def __str__(self):
        xlogical_str = ", ".join(str(ps) for ps in self.xlogical)
        zlogical_str = ", ".join(str(ps) for ps in self.zlogical)
        return f"X-Logical: [{xlogical_str}]\nZ-Logical: [{zlogical_str}]"

    def print_tableau_num(self):
        print("X-Logical:")
        for ps in self.xlogical:
            row = [f"{xpow} {zpow}" for xpow, zpow in zip(ps.xpow, ps.zpow)]
            print("\t".join(row))
        print("Z-Logical:")
        for ps in self.zlogical:
            row = [f"{xpow} {zpow}" for xpow, zpow in zip(ps.xpow, ps.zpow)]
            print("\t".join(row))

    def identity(self):
        """
        Creates a Tableau representing the identity operator.
        """
        for i in range(self.num_qudits):
            self.xlogical[i][i] = "X"
            self.zlogical[i][i] = "Z"

    def gate1(self, xmap, zmap):
        """
        Creates a Tableau representing a single qudit gate.
        Args:
            xmap: The output-side observable assuming the input-side is the logical X operator
            zmap: The output-side observable assuming the input-side is the logical Z operator
        """
        self.xlogical[0] = str(xmap)
        self.zlogical[0] = str(zmap)
        return self

    def gate2(self, xmap, zmap):
        """
        Creates a Tableau representing a two-qudit gate.
        Args:
            xmap: The output-side observable assuming the input-side is the logical X operator
            zmap: The output-side observable assuming the input-side is the logical Z operator
        """
        self.xlogical = [PauliString(self.num_qudits, x) for x in xmap]
        self.zlogical = [PauliString(self.num_qudits, z) for z in zmap]


class Program:
    def __init__(self, circuit, tableau=None):
        if tableau is None:
            self.stabilizer_tableau = Tableau(circuit.num_qudits, circuit.dimension)
            self.stabilizer_tableau.identity()
        else:
            self.stabilizer_tableau = tableau
        self.circuit = circuit

    def simulate(self):
        for time, gate in enumerate(self.circuit.operations):
            print("Time step", time)
            self.stabilizer_tableau.print_tableau_num()
            self.stabilizer_tableau = apply_gate(self.stabilizer_tableau, gate)
            print("\n")
        print("Final state")
        self.stabilizer_tableau.print_tableau_num()
