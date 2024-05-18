#ifndef QUANTUM_H
#define QUANTUM_H
#include <immintrin.h>
#include <vector>
struct QProg
{
	int num_qubits;         // # of qubits
	int gate_count;         // # of gates
	std::vector<char> op; // Instruction opcode
	std::vector<int> control_qubit; // Qubit 1
	std::vector<int> target_qubit; // Qubit 2 (target for CNOT)
	bool DISPQSTATE; // whether to print the state (q for final state only, Q for every iteration)
	bool DISPTIME; // whether to print the execution time
	bool SILENT;         // whether NOT to print measurement results
	bool DISPPROG; // whether to print instructions being executed as they're executed
	bool SUPPRESSM; // whether to suppress actual computation of determinate measurement results
};

struct QState
{
	int num_qubits;         // # of qubits
	std::vector<std::vector<__m512i>> x_bits; // (2n+1)*n matrix for stabilizer/destabilizer x bits (there's one "scratch row" at
	std::vector<std::vector<__m512i>> z_bits; // (2n+1)*n matrix for z bits                                                 the bottom)
	std::vector<int> phase;         // Phase bits: 0 for +1, 1 for i, 2 for -1, 3 for -i.  Normally either 0 or 2.
	unsigned int pw[32]; // pw[i] = 2^i
	int over512; // floor(n/512)+1
};
#endif // QUANTUM_H