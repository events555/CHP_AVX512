#ifndef GATES_H
#define GATES_H
void hadamard(struct QState *q, int b);
void phase(struct QState *q, long b);
void measure(struct QState *q, long b, int suppress);
void cnot(struct QState *q, long c, long t);
#endif