#ifndef __PDESolver1D__
#define __PDESolver1D__

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

class PDESolver1D{

    public:

    PDESolver1D(double, double, double, double(*func3)(double ), double(*func4)(double),void(*func1)(int, double*,double*), void(*func2)(int, double*,double*)); //constructor
    ~PDESolver1D(); //Destructor
    void Solve(double, double, bool); //solver (unparallelized)
	//void PointConvergenceTest(double (*func1)(double), double (*func2)(double), string);//Tests the point-wise convergence of the used Fin. Diff. Method
    //void NormConvergenceTest(double (*func1)(double), double (*func2)(double), string); //Tests norm convergence

    private:

    int N;
    double x_min;
    double x_max;
    double dx;
    double* Axis;
    double* y;
    void (*RHS1) (int, double*, double*);
    void (*RHS2) (int, double*, double*);
    double (*Initial_Condition_Phi) (double);
    double (*Initial_Condition_Pi) (double);

    void SetInitialConditions();
    void RHS(double*, double* );
    void RuggeKutta(double);
    void Saving_Data_File(int);
    //void Info_Conv_Test(double, string);

};

#endif
