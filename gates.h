#ifndef GATES_H
#define GATES_H
void hadamard(struct QState *q, long b);
void phase(struct QState *q, long b);
void measure(struct QState *q, long b, int suppress);
#endif