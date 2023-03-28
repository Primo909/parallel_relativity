//Advanced Topics of Computational Physics (2022/2023)
//Project 2 - The Parallelized Solution of Time Evolution PDEs
//by
//Catarina Corte-Real ist191035
//Francisco Valério Raposo ist196531
//Kevin Steiner ist1107611

//Wave equation

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
#include <string>
#include <chrono>

#include "header.h"
#include "WaveEquationSolver1D.h"
#include "NonLinearWaveEquationSolver1D.h"

using namespace std;
using namespace std::chrono;


double dx=1E-2; //Space discretization (uniform)
const double dt=1E-3; //Time discret.
const double x_min=-1, x_max=1; //Space interval 
const double simul_time = 2; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
int N=(x_max - x_min) / dx; //Number of spatial cells
const int number_ghosts = 3; // Number of ghost cells

double sigma = 0.1;
double x0 = 0;


double Sin(double x){
  return sin(2*M_PI*x);
}
double Zero(double x){
  return 0;
}
double GaussianFixed(double x){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}

int main(){

   // define fields
  double* y = new double[2*N];
  double* phi = &y[0];
  double* pi = &y[N];

  // impose initial values on phi and pi
  double x_aux; 
  
  // for a sinusoidal
  /*for(int j=0; j<N; j++){
    x_aux = x_min+j*dx;
    phi[j] = sin(2*M_PI*x_aux);
    pi[j] = sin(2*M_PI*x_aux);
  }*/
  // for a Gaussian
  double x0=0;
	double sigma=0.5;
  for(int j=0; j<N; j++){
    x_aux = x_min+j*dx;

    phi[j] = GaussianFixed(x_aux);
    pi[j] = 0;
  }
  WaveEquationSolver1D WaveEquation(-1,1,GaussianFixed,Zero);
  WaveEquation.Solve(dx,dt,simul_time, "a");

  /*NonLinearWaveEquationSolver1D NLWaveEquation;
  NLWaveEquation.Solve(x_min,x_max,dx,y,dt,simul_time,true);*/

  //WaveEquation.PointConvergenceTest(Sin, Zero, "./Data/sin_zero.dat");
  //WaveEquation.PointConvergenceTest(GaussianFixed, Zero, "./Data/gauss_zero.dat");
  //WaveEquation.NormConvergenceTest(GaussianFixed, Zero, "./Data/norm_conv_gauss_zero.dat");

  return 0;
  
}
