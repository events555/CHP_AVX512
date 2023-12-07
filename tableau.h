// functions.h

#ifndef TABLEAU_H
#define TABLEAU_H

void rowcopy(QState *q, long i, long k);
void rowswap(QState *q, long i, long k);
void rowset(QState *q, long i, long b);
void rowmult(QState *q, long i, long k);
int clifford(QState *q, long i, long k);

#endif
