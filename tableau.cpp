#include <immintrin.h>
#include <omp.h>
#include "quantum.cpp"
#include "tableau.h"

void rowcopy(QState *q, int i, int k) {

    for (int j = 0; j < q->over512; j++)
    {
        q->x_bits[i][j] = q->x_bits[k][j];
        q->z_bits[i][j] = q->z_bits[k][j];
    }
    q->phase[i] = q->phase[k];
    return;
}

void rowswap(QState *q, int i, int k)
{
    rowcopy(q, 2*q->num_qubits, k);
    rowcopy(q, k, i);
    rowcopy(q, i, 2*q->num_qubits);
    return;
}

void rowset(QState *q, int i, int b)
{

    unsigned int bit_array[16];

    for (int j = 0; j < q->over512; j++)
    {
        q->x_bits[i][j] = _mm512_setzero_si512();
        q->z_bits[i][j] = _mm512_setzero_si512();
    }
    q->phase[i] = 0;
    if (b < q->num_qubits)
    {
        int col_index = b >> 9;
        int int_index = (b & 511) >> 5;
        _mm512_storeu_si512(bit_array, q->x_bits[i][col_index]);         
        bit_array[int_index] = q->pw[b&31];
        q->x_bits[i][col_index] = _mm512_loadu_si512(bit_array);
    }
    else
    {
        int col_index = (b - q->num_qubits) >> 9;
        int int_index = ((b - q->num_qubits)&511) >> 5;
        _mm512_storeu_si512(bit_array, q->z_bits[i][col_index]);         
        bit_array[int_index] = q->pw[(b-q->num_qubits)&31];
        q->z_bits[i][col_index] = _mm512_loadu_si512(bit_array);
    }
    return;
}

int clifford(QState *q, int i, int k) {
    unsigned int pw;
    int e=0; // Power to which i is raised
    int i_int_index = (i & 511) >> 5;
    int k_int_index = (k & 511) >> 5;
    unsigned int k_x_bit_array[16];
    unsigned int k_z_bit_array[16];
    unsigned int i_x_bit_array[16];
    unsigned int i_z_bit_array[16];

    for (int j = 0; j < q->over512; j++)
    {
        for (int l = 0; l < 32; l++){
            pw = q->pw[l];
            _mm512_storeu_si512(k_x_bit_array, q->x_bits[k][j]);
            _mm512_storeu_si512(k_z_bit_array, q->z_bits[k][j]);
            _mm512_storeu_si512(i_x_bit_array, q->x_bits[i][j]);
            _mm512_storeu_si512(i_z_bit_array, q->z_bits[i][j]);
            if ((k_x_bit_array[k_int_index]&pw) && !(k_z_bit_array[k_int_index]&pw)){
                if((i_x_bit_array[i_int_index]&pw) && (i_z_bit_array[i_int_index]&pw)) e++;
                if(!(i_x_bit_array[i_int_index]&pw) && (i_z_bit_array[i_int_index]&pw)) e--;
            }
            if ((k_x_bit_array[k_int_index]&pw) && (k_z_bit_array[k_int_index]&pw)){
                if(!(i_x_bit_array[i_int_index]&pw) && (i_z_bit_array[i_int_index]&pw)) e++;
                if((i_x_bit_array[i_int_index]&pw) && !(i_z_bit_array[i_int_index]&pw)) e--;
            }
            if (!(k_x_bit_array[k_int_index]&pw) && (k_z_bit_array[k_int_index]&pw)){
                if((i_x_bit_array[i_int_index]&pw) && !(i_z_bit_array[i_int_index]&pw)) e++;
                if((i_x_bit_array[i_int_index]&pw) && (i_z_bit_array[i_int_index]&pw)) e--;
            }
        }
    }
    e = (e+q->phase[i]+q->phase[k])%4;
    if (e>=0) 
        return e;
    else 
        return e+4;
}

void rowmult(QState *q, int i, int k)
{

    q->phase[i] = clifford(q,i,k);
    for (int j = 0; j < q->over512; j++)
    {
        q->x_bits[i][j] = _mm512_xor_si512(q->x_bits[i][j], q->x_bits[k][j]);
        q->z_bits[i][j] = _mm512_xor_si512(q->z_bits[i][j], q->z_bits[k][j]);
    }
    return;
}
