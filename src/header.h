#ifndef __HEADER__
#define __HEADER__

void SecondDerivative(double*, double*);
void RHS(double*, double* );
void RungeKutta(double* , void (*Func)(double*, double*));

#endif