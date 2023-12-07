#include <omp.h>
#include "tableau.h"
#include "quantum.cpp"

void rowcopy(struct QState *q, long i, long k)
{
    long j;
    #pragma omp parallel for
    for (j = 0; j < q->over32; j++)
    {
        q->x[i][j] = q->x[k][j];
        q->z[i][j] = q->z[k][j];
    }
    q->r[i] = q->r[k];
    return;
}

void rowswap(struct QState *q, long i, long k)
{
    rowcopy(q, 2*q->n, k);
    rowcopy(q, k, i);
    rowcopy(q, i, 2*q->n);
    return;
}

void rowset(struct QState *q, long i, long b)
{
    long j;
    long b5;
    unsigned long b31;

    #pragma omp parallel for
    for (j = 0; j < q->over32; j++)
    {
        q->x[i][j] = 0;
        q->z[i][j] = 0;
    }
    q->r[i] = 0;
    if (b < q->n)
    {
        b5 = b>>5;
        b31 = b&31;
        q->x[i][b5] = q->pw[b31];
    }
    else
    {
        b5 = (b - q->n)>>5;
        b31 = (b - q->n)&31;
        q->z[i][b5] = q->pw[b31];
    }
    return;
}

int clifford(QState *q, long i, long k)
{
    long j;
    long l;
    unsigned long pw;
    long e=0; // Power to which i is raised

    #pragma omp parallel for reduction(+:e)
    for (j = 0; j < q->over32; j++)
    {
        for (l = 0; l < 32; l++)
        {
            pw = q->pw[l];
            if ((q->x[k][j]&pw) && (!(q->z[k][j]&pw))) // X
            {
                    if ((q->x[i][j]&pw) && (q->z[i][j]&pw)) e++;         // XY=iZ
                    if ((!(q->x[i][j]&pw)) && (q->z[i][j]&pw)) e--;         // XZ=-iY
            }
            if ((q->x[k][j]&pw) && (q->z[k][j]&pw))                                 // Y
            {
                    if ((!(q->x[i][j]&pw)) && (q->z[i][j]&pw)) e++;         // YZ=iX
                    if ((q->x[i][j]&pw) && (!(q->z[i][j]&pw))) e--;         // YX=-iZ
            }
            if ((!(q->x[k][j]&pw)) && (q->z[k][j]&pw))                         // Z
            {
                    if ((q->x[i][j]&pw) && (!(q->z[i][j]&pw))) e++;         // ZX=iY
                    if ((q->x[i][j]&pw) && (q->z[i][j]&pw)) e--;         // ZY=-iX
            }
        }
    }
    e = (e+q->r[i]+q->r[k])%4;
    if (e>=0) return e;
    else return e+4;
}

void rowmult(struct QState *q, long i, long k)
{
    long j;

    q->r[i] = clifford(q,i,k);
    #pragma omp parallel for
    for (j = 0; j < q->over32; j++)
    {
        q->x[i][j] ^= q->x[k][j];
        q->z[i][j] ^= q->z[k][j];
    }
    return;
}
