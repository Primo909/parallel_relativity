#ifndef __WaveEquationSolver1D__
#define __WaveEquationSolver1D__

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
#include <string>
#include <chrono>
//#include "mpi.h"

#include "header.h"

using namespace std;
using namespace std::chrono;

double Gaussian(double, double, double);

class WaveEquationSolver1D{

    public:

    WaveEquationSolver1D(double, double, double(*func1)(double ), double(*func2)(double)); //constructor
    void Solve(double, double, double, string); //solver (unparallelized)
    //void Parallel_Solve(double, double, double, string); //solver (unparallelized)
	void PointConvergenceTest(string);//Tests the point-wise convergence of the used Fin. Diff. Method
    void NormConvergenceTest(string); //Tests norm convergence

    private:

    double x_max;
    double x_min;
    double (*Initial_Condition_Phi) (double);
    double (*Initial_Condition_Pi) (double);
    //int size;
    //int rank;

    //
    void SetInitialConditions(int, double*, double*);
    void FirstDerivative();
    void SecondDerivative(int, double, double* , double* );
    void RHS(int, double, double*, double* );
    void RungeKutta(int, double, double, double*);
    void Saving_Data_File(int ,double* , double*, int );
    void Info_Conv_Test(double, string);

    /*//Parallel Computation
    void Parallel_SetInitialConditions(int, double*, double*);
    void Parallel_RungeKutta(int, double, double, double*);
    void Parallel_RHS(int, double, double*, double* );
    void Parallel_SecondDerivative(int, double, double* , double* );
    */

};

#endif
