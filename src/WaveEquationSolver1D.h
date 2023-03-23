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

using namespace std;
using namespace std::chrono;

const string Simulation_Start_String= "Simulation has started";
const string Simulation_End_String= "Simulation has ended";
const string File_String = "./Data/numerical_phi_ite";

double Gaussian(double, double, double);

class WaveEquationSolver1D{

    public:

    WaveEquationSolver1D(); //constructor
    void Solve(double, double, double, double*, double, double); //solver
	void PointConvergenceTest(double (*func1)(double), double (*func2)(double), string);//Tests the convergence of the used Fin. Diff. Method

    private:

    void SetInitialConditions(int, double*, double*);
    void FirstDerivative();
    void SecondDerivative(int, double, double* , double* );
    void RHS(int, double, double*, double* );
    void RuggeKutta(int, double, double, double*);
    void Saving_Data_File(int ,double* , double*, int );

};

#endif
