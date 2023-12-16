#include <immintrin.h>
#include <cstdint>
#include "quantum.cpp"
void hadamard(QState *q, int b)
{
    __m512i tmp;
    int b6 = b / 64;
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

void phase(QState *q, int b)
{
    int b6 = b / 64;
    uint64_t pw = q->pw[b % 64];
    __m512i pw_vec = _mm512_set1_epi64(pw);

    for (int i = 0; i < 2 * q->num_qubits; i++)
    {
        if (_mm512_test_epi64_mask(q->x_bits[i][b6], pw_vec) && _mm512_test_epi64_mask(q->z_bits[i][b6], pw_vec))
            q->phase[i] = (q->phase[i] + 2) % 4;

        q->z_bits[i][b6] = _mm512_xor_si512(q->z_bits[i][b6], _mm512_and_si512(q->x_bits[i][b6], pw_vec));
    }
}

void cnot(struct QState *q, int b, int c)
{
    int control_index = b >> 9;
    int target_index = c >> 9;
    uint64_t pwb = q->pw[b & 31];
    uint64_t pwc = q->pw[c & 31];

    for (int i = 0; i < 2 * q->num_qubits; i++)
    {
        uint64_t val_x_control[8];
        uint64_t val_z_control[8];
        uint64_t val_x_target[8];
        uint64_t val_z_target[8];
        _mm512_storeu_si512(val_x_control, q->x_bits[i][control_index]);
        _mm512_storeu_si512(val_z_control, q->z_bits[i][control_index]);
        _mm512_storeu_si512(val_x_target, q->x_bits[i][target_index]);
        _mm512_storeu_si512(val_z_target, q->z_bits[i][target_index]);

        //To do: extract the bit from the array 

        // _mm512_store_si512((__m512i *)q->x_bits[i][target_index], _mm512_xor_si512(target_x_bits_xor_pwc, control_x_bits_and_pwb));
        // _mm512_store_si512((__m512i *)q->z_bits[i][control_index], _mm512_xor_si512(control_z_bits_and_pwb, target_z_bits_xor_pwc));

        // __mmask8 x_bits_mask = _mm512_test_epi64_mask(control_x_bits, pwb_vec);
        // __mmask8 z_bits_mask = _mm512_test_epi64_mask(target_z_bits, pwc_vec);

        // if (x_bits_mask && z_bits_mask)
        //     q->phase[i] = (q->phase[i] + 2) % 4;

        // if (x_bits_mask && z_bits_mask && _mm512_testz_epi64_mask(target_x_bits, pwc_vec) && _mm512_testz_epi64_mask(control_z_bits, pwb_vec))
        //     q->phase[i] = (q->phase[i] + 2) % 4;
    }
}