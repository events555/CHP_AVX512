#include <immintrin.h>
#include <cstdint>
#include "quantum.cpp"
void hadamard(QState *q, long b)
{
    __m512i tmp;
    long b6 = b / 64;
    uint64_t pw = q->pw[b % 64];
    __m512i pw_vec = _mm512_set1_epi64(pw);

    for (int i = 0; i < 2 * q->num_qubits; i++)
    {
        tmp = q->x_bits[i][b6];
        q->x_bits[i][b6] = _mm512_xor_si512(q->x_bits[i][b6], _mm512_and_si512(_mm512_xor_si512(q->x_bits[i][b6], q->z_bits[i][b6]), pw_vec));
        q->z_bits[i][b6] = _mm512_xor_si512(q->z_bits[i][b6], _mm512_and_si512(_mm512_xor_si512(q->z_bits[i][b6], tmp), pw_vec));

        if (_mm512_test_epi64_mask(q->x_bits[i][b6], pw_vec) && _mm512_test_epi64_mask(q->z_bits[i][b6], pw_vec))
            q->phase[i] = (q->phase[i] + 2) % 4;
    }
}

void phase(QState *q, long b)
{
    long b6 = b / 64;
    uint64_t pw = q->pw[b % 64];
    __m512i pw_vec = _mm512_set1_epi64(pw);

    for (int i = 0; i < 2 * q->num_qubits; i++)
    {
        if (_mm512_test_epi64_mask(q->x_bits[i][b6], pw_vec) && _mm512_test_epi64_mask(q->z_bits[i][b6], pw_vec))
            q->phase[i] = (q->phase[i] + 2) % 4;

        q->z_bits[i][b6] = _mm512_xor_si512(q->z_bits[i][b6], _mm512_and_si512(q->x_bits[i][b6], pw_vec));
    }
}
