// functions.h

#ifndef TABLEAU_H
#define TABLEAU_H

void rowcopy(QState *q, int i, int k);
void rowswap(QState *q, int i, int k);
void rowset(QState *q, int i, int b);
void rowmult(QState *q, int i, int k);
int clifford(QState *q, int i, int k);

#endif
