import numpy as np
import mats
import os
import math

SUPPRESS_OUTPUT = False
RUN_ROUNDS_FROM_1_TO_N = True
LOG_OUTPUT = True
imag_unit = complex(0+1j)

def qudit_stab(qudits, gates, rounds):
    #initialize dict
    vectors = {}
    for q in qudits:
        vectors[q.tobytes()] = np.copy(q)
    
    for r in range(rounds):
        for p in gates:
            curr_vectors = list(vectors.values())
            for v in curr_vectors:
                vec = np.matmul(p, v)
                vectors[vec.tobytes()] = np.copy(vec)

    return vectors

#much slower than using dict, but you can use allclose for matrix comparison
def qudit_stab_list(qudits, gate_seq, gates, symb, rounds):

    vectors = []
    d = qudits[0].shape[0]

    for q in qudits:
        vectors.append(q)

    for r in range(rounds):
        curr_vectors = vectors.copy()
        for vec_ind, v in enumerate(curr_vectors):
             for gate_ind, g in enumerate(gates):
                 vec = np.matmul(g, v)
                 not_close = True
                 omega = np.exp((-2 * math.pi * imag_unit) / (2 * d))
                 for k in range(2 * d):
                     if any( [ np.allclose(u, np.power(omega, k) * vec) for u in vectors ] ):
                         not_close = False
                         break

                 if not_close:
                     vectors.append(vec)
                     gate_seq.append(symb[gate_ind] + gate_seq[vec_ind]) 
    
    return [vectors, gate_seq]
       


def qudit_comp_basis(d):
    indicator = lambda n, k: [1+0j] if k == n else [0+0j]
    return np.array([ [indicator(dd, ddd) for ddd in range(d)] for dd in range(d)])


if __name__ == "__main__":
    
    #sanitize input later
    d = int(input("Dim: "))
    
    #gates = mats.pauli_matrices if d == 2 else mats.chp_3[:7]
    
    #initialization 
    gates, symb = mats.chp_plus_labels(d)
    #ignore CNOT/SUM for now
    symb = symb[1:]
    gates = gates[1:]
    gate_seq = ["|" + str(k) + ">" for k in range(d)]
    prev_gates = 0
    stab = qudit_comp_basis(d)
    rounds_needed = 0
    #initialization

    while True:
        stab, gate_seq = qudit_stab_list(stab, gate_seq, gates, symb, 1)
        rounds_needed += 1

        if not SUPPRESS_OUTPUT:
            for v in stab:
                print(v)
                print()
                
        num_gates = len(stab)
        
        if (prev_gates == num_gates):

            print("Total gates seen: " + str(num_gates))
            if LOG_OUTPUT:
                fname = "d_" + str(d) + "_stabilized_states.txt"

                if os.path.exists(fname):
                    os.remove(fname)

                fp = open(fname, "w")
                
                for v_ind, v in enumerate(stab):
                    fp.write(str(gate_seq[v_ind]) + "\n")
                    fp.write(str(v) + "\n\n")

                print("Wrote gates to " + str(fname))
                fp.close()

            break

        prev_gates = num_gates
        print("The number of seen vectors: " + str(len(stab)))
        print()

    print("Successfully terminated in " + str(rounds_needed) + " rounds.")
