#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>


#include "Grid1D.h"

#ifndef __Var1D__
#define __Var1D__

using namespace std;

class Var1D{
    
    public:

    Var1D(GRID1D, float*); //constructor - receives array with spatial values
    void SetGhostCells(); //sets ghost cells (with periodic boundary conditions)
    void Derivative(float *, int ); //calculates 1st derivative (order O)
    void SecondDerivative(float *, int ); //calculates second derivative (order O)

    private:

    int N;
    float dx;
    float* V; //array with spatial values
    float* dVdx; //first derivative
    float* ddVddx; //second derivative
    float* RightGhostCells; //array with ghost cells on the right
    float* LeftGhostCells; //array with ghost cells in the left
    int NGhostCells=3; //# of ghost cells (same for either side)

};

#endif