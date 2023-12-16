#include <immintrin.h>
#include <omp.h>
#include "quantum.cpp"
#include "tableau.h"

void rowcopy(QState *q, int i, int k) {

    #pragma omp parallel for
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

    int bit_array[16];

    for (int j = 0; j < q->over512; j++)
    {
        q->x_bits[i][j] = _mm512_set1_epi64(0);
        q->z_bits[i][j] = _mm512_set1_epi64(0);
    }
    q->phase[i] = 0;
    if (b < q->num_qubits)
    {
        int col_index = b >> 9;
        int int_index = b >> 5;
        _mm512_storeu_si512(bit_array, q->x_bits[i][col_index]);         
        bit_array[int_index] = q->pw[int_index];
        q->x_bits[i][col_index] = _mm512_loadu_si512(bit_array);
    }
    else
    {
        int col_index = (b - q->num_qubits) >> 9;
        int int_index = (b - q->num_qubits) >> 5;
        _mm512_storeu_si512(bit_array, q->z_bits[i][col_index]);         
        bit_array[int_index] = q->pw[int_index];
        q->z_bits[i][col_index] = _mm512_loadu_si512(bit_array);
    }
    return;
}

int clifford(QState *q, int i, int k)
{
    unsigned int pw;
    int e=0; // Power to which i is raised

    #pragma omp parallel for reduction(+:e)
    for (int j = 0; j < q->over512; j++)
    {
        for (int l = 0; l < 32; l++)
        {
            unsigned int array[32] = {0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0};
            array[l] = q->pw[l];
            __m512i pwmask = _mm512_loadu_si512(array);
            
            __m512i x_and_pw_k = _mm512_and_si512(q->x_bits[k][j], pwmask);
            __m512i z_and_pw_k = _mm512_and_si512(q->z_bits[k][j], pwmask);
            __m512i x_and_pw_i = _mm512_and_si512(q->x_bits[i][j], pwmask);
            __m512i z_and_pw_i = _mm512_and_si512(q->z_bits[i][j], pwmask);

            if (_mm512_test_epi32_mask(x_and_pw_k, _mm512_set1_epi32(-1)) && !_mm512_test_epi32_mask(z_and_pw_k, _mm512_set1_epi32(-1))) // x_bits
            {
                if (_mm512_test_epi32_mask(x_and_pw_i, _mm512_set1_epi32(-1)) && _mm512_test_epi32_mask(z_and_pw_i, _mm512_set1_epi32(-1))) e++;         // XY=iZ
                if (!_mm512_test_epi32_mask(x_and_pw_i, _mm512_set1_epi32(-1)) && _mm512_test_epi32_mask(z_and_pw_i, _mm512_set1_epi32(-1))) e--;         // XZ=-iY
            }
            if (_mm512_test_epi32_mask(x_and_pw_k, _mm512_set1_epi32(-1)) && _mm512_test_epi32_mask(z_and_pw_k, _mm512_set1_epi32(-1)))                                 // Y
            {
                if (!_mm512_test_epi32_mask(x_and_pw_i, _mm512_set1_epi32(-1)) && _mm512_test_epi32_mask(z_and_pw_i, _mm512_set1_epi32(-1))) e++;         // YZ=iX
                if (_mm512_test_epi32_mask(x_and_pw_i, _mm512_set1_epi32(-1)) && !_mm512_test_epi32_mask(z_and_pw_i, _mm512_set1_epi32(-1))) e--;         // YX=-iZ
            }
            if (!_mm512_test_epi32_mask(x_and_pw_k, _mm512_set1_epi32(-1)) && _mm512_test_epi32_mask(z_and_pw_k, _mm512_set1_epi32(-1)))                         // z_bits
            {
                if (_mm512_test_epi32_mask(x_and_pw_i, _mm512_set1_epi32(-1)) && !_mm512_test_epi32_mask(z_and_pw_i, _mm512_set1_epi32(-1))) e++;         // ZX=iY
                if (_mm512_test_epi32_mask(x_and_pw_i, _mm512_set1_epi32(-1)) && _mm512_test_epi32_mask(z_and_pw_i, _mm512_set1_epi32(-1))) e--;         // ZY=-iX
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
    #pragma omp parallel for
    for (int j = 0; j < q->over512; j++)
    {
        q->x_bits[i][j] ^= q->x_bits[k][j];
        q->z_bits[i][j] ^= q->z_bits[k][j];
    }
    return;
}
