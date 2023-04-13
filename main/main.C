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

double dx = 0;//Space discretization (uniform)
const double dt=1E-4; //Time discret.
const double x_min=-1, x_max=1; //Space interval 
const double simul_time = 2; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
int N=1000;  //Number of spatial cells
const int number_ghosts = 2; // Number of ghost cells

double x0 = 0;
double sigma = 0.5;  

double GaussianFixed(double x){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}

double Sin(double x){
  return sin(2*M_PI*x);
}
double Zero(double x){
  return 0;
}


int main(int argc, char* argv[]){

  // DELETE THIS LATER ON 
  int N = atoi(argv[1]); 
  double dx=(x_max-x_min)/N;
  cout << "dx = " << dx << endl;

  // counting time
  auto start = high_resolution_clock::now(); 

  // define fields
  double* y = new double[2*N];
  double* phi = &y[0];
  double* pi = &y[N];

  // solving a wave equation

  WaveEquationSolver1D WaveEquation(-1,1,GaussianFixed,Zero);

  WaveEquation.Solve(dx,dt,simul_time, "");

  // stop counting time
  auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
  cout<<Simulation_End_String<<endl;

  // DELETE THIS LATER ON
  cout << "control," << N << "," << duration.count()/1E6 << endl;
  
  return 0;
  
}
