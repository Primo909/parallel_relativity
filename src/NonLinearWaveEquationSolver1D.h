#ifndef __NonLinearWaveEquationSolver1D__
#define __NonLinearWaveEquationSolver1D__

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
#include <string>
#include <chrono>

#include "header.h"

using namespace std;
using namespace std::chrono;

double Gaussian(double, double, double);

class NonLinearWaveEquationSolver1D{

    public:

    NonLinearWaveEquationSolver1D(); //constructor
    void Solve(double, double, double, double*, double, double, bool); //solver (unparallelized)
	void PointConvergenceTest(double (*func1)(double), double (*func2)(double), string);//Tests the point-wise convergence of the used Fin. Diff. Method
    void NormConvergenceTest(double (*func1)(double), double (*func2)(double), string); //Tests norm convergence

    private:

    void SetInitialConditions(int, double*, double*);
    void FirstDerivative();
    void SecondDerivative(int, double, double* , double* );
    void RHS(int, double, double*, double* );
    void RuggeKutta(int, double, double, double*);
    void Saving_Data_File(int ,double* , double*, int );
    void Info_Conv_Test(double, string);

};

#endif