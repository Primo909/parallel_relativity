//Advanced Topics of Computational Physics (2022/2023)
//Project 2 - The Parallelized Solution of Time Evolution PDEs
//by
//Catarina Corte-Real ist191035
//Francisco Valério Raposo ist19651
//Kevin Steiner ist1107611

//One dimensional standard wave equation
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 

#include "header.h"

using namespace std;

const double dx=1E-3; //Space discretization (uniform)
const double dt=1E-3; //Time discret.
const double x_min=0, x_max=1; //Space interval 
const double simul_time = 0.1; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
const int N=(x_max - x_min) / dx; //Number of spatial cells
const int number_ghosts = 3; // Number of ghost cells


void SecondDerivative(double* field, double* second_derivative_vector){

  // creates ghost cells
  double* right_ghost_cells_field = new double[number_ghosts];
  double* left_ghost_cells_field = new double[number_ghosts];
  // assign values to ghost cells
  for(int j=0; j<number_ghosts; j++){
    right_ghost_cells_field[j] = field[j];
    left_ghost_cells_field[number_ghosts-1-j] = field[N-1-j];
  }

  // calculate the second derivative with five points (-2,-1,0,1,2)
  // formula from https://web.media.mit.edu/~crtaylor/calculator.html

  for(int j=0; j<N; j++){
    if(j==0) second_derivative_vector[j] = (-1*left_ghost_cells_field[number_ghosts-2]+16*left_ghost_cells_field[number_ghosts-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
    else if(j==1) second_derivative_vector[j] = (-1*left_ghost_cells_field[number_ghosts+1-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
    else if(j==N-1) second_derivative_vector[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*right_ghost_cells_field[0]-1*right_ghost_cells_field[1])/(12*1.0*dx*dx);
    else if(j==N-2) second_derivative_vector[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*right_ghost_cells_field[0])/(12*1.0*dx*dx);
    else second_derivative_vector[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
  }
  // delete the ghost cells
  delete[] left_ghost_cells_field, right_ghost_cells_field;
}

//Right Hand Side
void RHS(double* state_vector, double* rhs_vector){
  double* phi = &state_vector[0];
  double* pi = &state_vector[N];
  for(int j=0; j<N; j++) rhs_vector[j] = pi[j];
  SecondDerivative(phi, &rhs_vector[N]);
  
  delete[] phi, pi;
}

//4th order Runge-Kutta
void RungeKutta(double* state_vector, void (*Func)(double*, double*)){
  
  double* k1 = new double[2*N];
  double* k2 = new double[2*N];
  double* k3 = new double[2*N];
  double* k4 = new double[2*N];

  double* vector_temp = new double[2*N];

  RHS(state_vector, k1); // Calc k1
  for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dx*k1[j]/2;
  RHS(vector_temp, k2); // Calc k2
  for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dx*k2[j]/2; 
  RHS(vector_temp, k3); // Calc k3
  for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dx*k3[j];
  RHS(vector_temp, k4);// Calc k4
  
  for(int j=0; j<2*N; j++) state_vector[j] = state_vector[j]+(k1[j] + 2*k2[j] + 2*k3[j] +k4[j]) / 6; // Overwrite the state vector

}

int main(){

  // define fields
  double* y = new double[2*N];
  double* phi = &y[0];
  double* pi = &y[N];

  // impose initial values on phi and pi
  double x_aux;
  for(int j=0; j<N; j++){
    x_aux = j*dx;
    phi[j] = sin(2*M_PI*x_aux);
    pi[j] = 2*M_PI*cos(2*M_PI*x_aux);
  }

  RungeKutta(y, RHS);

  /* check if derivative works
  
  double* second_phi = new double[N];
  SecondDerivative(phi, second_phi);
  double* rhs_vector = new double[2*N];
  rhs(y, rhs_vector);

  cout << setprecision(4) << "j" << setw(12) << "phi" << setw(12) << "phi_xx" << endl;
  for(int j=0;j<N;j++){
  cout << setprecision(4) << j << setw(12) << phi[j] << setw(12) << -1*second_phi[j]/(4*M_PI*M_PI) << setw(12) << phi[j]+1*second_phi[j]/(4*M_PI*M_PI)<< endl;
  }
  */


  /*
  // check for f flip
  cout << setprecision(2) << "j" << setw(12) << "y" << setw(12) << "f" << endl;
  for(int j=0;j<2*N;j++){
  cout << setprecision(2) << j << setw(12) << y[j] << setw(12) << rhs_vector[j] << endl;
  }
  */
  // Runge kutta
  /*
  double* k1 = new double[2*N];
  double* k2 = new double[2*N];
  double* k3 = new double[2*N];
  double* k4 = new double[2*N];
  */
  
  // output something
  //cout << "j" << "," << "phi" << "," << "pi" << endl;
  //for(int j=0;j<N;j++){
  //cout << j << "," << phi[j] << "," << pi[j] << endl;
  //}
  // delete pointer arrays
  delete[] phi, pi, y;
    
  return 0;
  
}
