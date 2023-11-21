class Qudit:
    def __init__(self, dimension):
        self.dimension = dimension

class QuditRegister:
    def __init__(self, name, dimension, size):
        self.name = name
        self.dimension = dimension
        self.size = size
        self.qudits = [Qudit(dimension) for _ in range(size)]

    def __str__(self):
        return f"{self.name}: Qudit Register of {self.size} qudits, each of dimension {self.dimension}"

# Example usage:
qr = QuditRegister("qr", 3, 5)
gamma = Qudit(3)
