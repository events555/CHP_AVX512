#ifndef GATES_H
#define GATES_H
void hadamard(struct QState *q, int b);
void phase(struct QState *q, int b);
void measure(struct QState *q, int b, int suppress);
void cnot(struct QState *q, int c, int t);
#endif