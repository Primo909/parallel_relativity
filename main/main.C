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

#include "mpi.h"
#include "header.h"
#include "WaveEquationSolver1D.h"
#include "NonLinearWaveEquationSolver1D.h"

using namespace std;
using namespace std::chrono;

double dx=1E-3; //Space discretization (uniform)
const double dt=1E-4; //Time discret.
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

int main(int argc, char* argv[]){

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

  //double* y= new double[N];

  WaveEquationSolver1D WaveEquation(-1,1,GaussianFixed,Zero);

  //MPI_Init(&argc, &argv);

  /*int size,rank;
  int my_first, my_last;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  my_first = rank* N/size;
  my_last = (rank+1)*N/size;

  cout<<rank<<" "<<size<<endl;

  cout<<my_first<<" "<<my_last<<endl;

  for(int j=my_first; j<my_last; j++){
    x_aux = x_min+j*dx;

    phi[j] = GaussianFixed(x_aux);
    pi[j] = 0;
  }*/
  
  WaveEquation.Solve(dx,dt,simul_time, "");
  //WaveEquation.Solve(dx,dt,simul_time, "");
  //WaveEquation.Solve(dx,dt,simul_time, "a");
  //WaveEquation.Parallel_Solve(dx,dt,simul_time, "");
  
  //WaveEquation.PointConvergenceTest("./Data/gauss_zero.dat");
  //WaveEquation.NormConvergenceTest("./Data/norm_conv_gauss_zero.dat");

  //MPI_Finalize();

  


  
  return 0;
  
}
