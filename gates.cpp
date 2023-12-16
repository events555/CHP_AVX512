#include <immintrin.h>
#include <cstdint>
#include "quantum.cpp"
#include "tableau.cpp"

void hadamard(QState *q, int b)
{
    int col_index = b >> 9;
    int int_index = b >> 5;
    unsigned int pw = q->pw[b&31];
    int x_bits_array[16];
    int z_bits_array[16];
    for (int i = 0; i < 2 * q->num_qubits; i++)
    {
        _mm512_storeu_si512(x_bits_array, q->x_bits[i][col_index]);
        _mm512_storeu_si512(z_bits_array, q->z_bits[i][col_index]);

         int tmp = x_bits_array[int_index];
         x_bits_array[int_index] ^= (x_bits_array[int_index] ^ z_bits_array[int_index]) & pw;
         z_bits_array[int_index] ^= (z_bits_array[int_index] ^ tmp) & pw;
         if ((x_bits_array[int_index]&pw) && (z_bits_array[int_index]&pw))
            q->phase[i] = (q->phase[i]+2)%4;

        q->x_bits[i][col_index] = _mm512_loadu_si512(x_bits_array);
        q->z_bits[i][col_index] = _mm512_loadu_si512(z_bits_array);
    }
}

void phase(QState *q, int b)
{
    int col_index = b >> 9;
    int int_index = b >> 5;
    unsigned int pw = q->pw[b&31];
    int x_bits_array[16];
    int z_bits_array[16];
    for (int i = 0; i < 2 * q->num_qubits; i++)
    {
        _mm512_storeu_si512(x_bits_array, q->x_bits[i][col_index]);
        _mm512_storeu_si512(z_bits_array, q->z_bits[i][col_index]);

        if ((x_bits_array[int_index]&pw) && (z_bits_array[int_index]&pw)) 
            q->phase[i] = (q->phase[i]+2)%4;
        z_bits_array[int_index] ^= x_bits_array[int_index]&pw;

        q->x_bits[i][col_index] = _mm512_loadu_si512(x_bits_array);
        q->z_bits[i][col_index] = _mm512_loadu_si512(z_bits_array);
    }
}

void cnot(struct QState *q, int b, int c)
{
    int control_index = b >> 9;
    int control_bit = b >> 5;
    int target_index = c >> 9;
    int target_bit = c >> 5;
    int control_x_bits_array[16];
    int control_z_bits_array[16];
    int target_x_bits_array[16];
    int target_z_bits_array[16];
    unsigned int pwb = q->pw[b&31];
    unsigned int pwc = q->pw[c&31];
	for (int i = 0; i < 2*q->num_qubits; i++)
	{
        _mm512_storeu_si512(control_x_bits_array, q->x_bits[i][control_index]);
        _mm512_storeu_si512(control_z_bits_array, q->z_bits[i][control_index]);
        _mm512_storeu_si512(target_x_bits_array, q->x_bits[i][target_index]);
        _mm512_storeu_si512(target_z_bits_array, q->z_bits[i][target_index]);

        if (control_x_bits_array[control_bit]&pwb) 
            target_x_bits_array[target_bit] ^= pwc;
        control_x_bits_array[control_bit] = target_x_bits_array[target_bit];
        if (target_z_bits_array[target_bit]&pwc) 
            control_z_bits_array[control_bit] ^= pwb; //RACE CONDITION: Control_z_bits is not the same as target_z_bits even though they may refer to same AVX register
        target_z_bits_array[target_bit] = control_z_bits_array[control_bit];
		if ((control_x_bits_array[control_bit]&pwb) && (target_z_bits_array[target_bit]&pwc) &&
			(target_x_bits_array[target_bit]&pwc) && (control_z_bits_array[control_bit]&pwb))
			    q->phase[i] = (q->phase[i]+2)%4;
		if ((control_x_bits_array[control_bit]&pwb) && (target_z_bits_array[target_bit]&pwc) &&
			!(target_x_bits_array[target_bit]&pwc) && !(control_z_bits_array[control_bit]&pwb))
				q->phase[i] = (q->phase[i]+2)%4;

        q->x_bits[i][control_index] = _mm512_loadu_si512(control_x_bits_array);
        q->z_bits[i][control_index] = _mm512_loadu_si512(control_z_bits_array);
        q->x_bits[i][target_index] = _mm512_loadu_si512(target_x_bits_array);
        q->z_bits[i][target_index] = _mm512_loadu_si512(target_z_bits_array);
	}
}

int measure(struct QState *q, int b, int sup)

// Measure qubit b
// Return 0 if outcome would always be 0
//                 1 if outcome would always be 1
//                 2 if outcome was random and 0 was chosen
//                 3 if outcome was random and 1 was chosen
// sup: 1 if determinate measurement results should be suppressed, 0 otherwise

{

	bool ran = 0;
	int i;
	int p; // pivot row in stabilizer
	int m; // pivot row in destabilizer
	int col_index = b>>9;
    int int_index = b>>5;
	unsigned int pw = q->pw[b&31];
    int x_bits_array[16];
	for (p = 0; p < q->num_qubits; p++)         // loop over stabilizer generators
	{
        _mm512_storeu_si512(x_bits_array, q->x_bits[p+q->num_qubits][col_index]);
        if (x_bits_array[int_index]&pw) 
            ran = true;         // if a Zbar does NOT commute with Z_b (the
        if (ran) break;                                                 // operator being measured), then outcome is random
	}

	// If outcome is indeterminate
	if (ran)
	{
        rowcopy(q, p, p + q->num_qubits);                         // Set Xbar_p := Zbar_p
        rowset(q, p + q->num_qubits, b + q->num_qubits);                 // Set Zbar_p := Z_b
        q->phase[p + q->num_qubits] = rand() % 2;                 // moment of quantum randomness
        for (i = 0; i < 2*q->num_qubits; i++) {                 // Now update the Xbar's and Zbar's that don't commute with
            _mm512_storeu_si512(x_bits_array, q->x_bits[i][col_index]);
                if ((i!=p) && (x_bits_array[int_index]&pw))         // Z_b
                        rowmult(q, i, p);
        }
        if (q->phase[p + q->num_qubits]) 
            return 3;
        else 
            return 2;
	}

	// If outcome is determinate
	if ((!ran) && (!sup))
	{
        for (m = 0; m < q->num_qubits; m++) {                        // Before we were checking if stabilizer generators commute
            _mm512_storeu_si512(x_bits_array, q->x_bits[m][col_index]);         // with Z_b; now we're checking destabilizer generators
            if (x_bits_array[int_index]&pw) break;                 // with Z_b; now we're checking destabilizer generators
        }
        rowcopy(q, 2*q->num_qubits, m + q->num_qubits);
        for (i = m+1; i < q->num_qubits; i++) {
            _mm512_storeu_si512(x_bits_array, q->x_bits[i][col_index]);
                if (x_bits_array[int_index]&pw)
                        rowmult(q, 2*q->num_qubits, i + q->num_qubits);
        }
        if (q->phase[2*q->num_qubits]) 
            return 1;
        else 
            return 0;
	}

	return 0;

}