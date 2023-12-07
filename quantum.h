#ifndef QUANTUM_H
#define QUANTUM_H
struct QProg
{
	long n;         // # of qubits
	long T;         // # of gates
	char *a; // Instruction opcode
	long *b; // Qubit 1
	long *c; // Qubit 2 (target for CNOT)
	int DISPQSTATE; // whether to print the state (q for final state only, Q for every iteration)
	int DISPTIME; // whether to print the execution time
	int SILENT;         // whether NOT to print measurement results
	int	DISPPROG; // whether to print instructions being executed as they're executed
	int SUPPRESSM; // whether to suppress actual computation of determinate measurement results
};

struct QState
{
	long n;         // # of qubits
	unsigned long **x; // (2n+1)*n matrix for stabilizer/destabilizer x bits (there's one "scratch row" at
	unsigned long **z; // (2n+1)*n matrix for z bits                                                 the bottom)
	int *r;         // Phase bits: 0 for +1, 1 for i, 2 for -1, 3 for -i.  Normally either 0 or 2.
	unsigned long pw[32]; // pw[i] = 2^i
	long over32; // floor(n/8)+1
};
#endif // QUANTUM_H