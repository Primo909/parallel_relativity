//Advanced Topics of Computational Physics (2022/2023)
//Project 2 - The Parallelized Solution of Time Evolution PDEs
//by
//Catarina Corte-Real ist191035
//Francisco Valério Raposo ist196531
//Kevin Steiner ist1107611

//One dimensional standard wave equation

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
#include <string>

#include "header.h"

using namespace std;

const double dx=1E-3; //Space discretization (uniform)
const double dt=1E-4; //Time discret.
const double x_min=0, x_max=1; //Space interval 
const double simul_time = 1; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
const int N=(x_max - x_min) / dx; //Number of spatial cells
const int number_ghosts = 3; // Number of ghost cells

const string File_String = "./Data/numerical_phi_ite";


void Saving_Data_File(double* state_vector, int iteration){

    string string_aux=File_String+to_string(iteration+1)+".dat";
    
    fstream ITERATION_FILE;
    ITERATION_FILE.open(string_aux,ios::out);

    if (!ITERATION_FILE){                 
        cout<<"Error: File not created"<<endl;    
    }else{
        cout<<"File created"<<endl;              
    }

    for(int j=0; j<N; j++) ITERATION_FILE<<j*dx<<" "<<state_vector[j]<<endl;

    ITERATION_FILE.close();   
}

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
}

//4th order Runge-Kutta
void RungeKutta(double* state_vector, void (*Func)(double*, double*)){

  double* k1 = new double[2*N];
  double* k2 = new double[2*N];
  double* k3 = new double[2*N];
  double* k4 = new double[2*N];

  double* vector_temp = new double[2*N];

  //cout<<"a"<<endl;

  Func(state_vector, k1); // Calc k1
  for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dt*k1[j]/2;
  Func(vector_temp, k2); // Calc k2
  for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dt*k2[j]/2; 
  Func(vector_temp, k3); // Calc k3
  for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dt*k3[j];
  Func(vector_temp, k4);// Calc k4
  for(int j=0; j<2*N; j++) state_vector[j] = state_vector[j]+(k1[j] + 2*k2[j] + 2*k3[j] +k4[j])*dt/6; // Overwrite the state vector

  delete[] vector_temp,k1,k2,k3,k4;
}

double Gaussian(double x, double sigma, double x0){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}

double Deviation(double result, double exact){
	return (result-exact)/exact;
}

void Convergence_Plot(double* state_vector, double* exact_phi, double dx_min, double dx_max, double num_it)
{
  ofstream CONVERGENCE_FILE("conv_file.dat");

    if (!CONVERGENCE_FILE){                 
        cout << "Error: File not created" << endl;    
    }else{
        cout << "File created" << endl;              
    }

  double* phi = &state_vector[0];
  double* pi = &state_vector[N];
  double step = (dx_max - dx_min)/num_it;
  double dx_i;

  for(int idx=0; idx < num_it; idx++)
  {
    dx_i = dx_min + idx*step;
    double x_aux; 
    double l1_norm = 0;

    for(int j=0; j<N; j++){
      x_aux = j*dx_i;
      // for a sinusoidal
      phi[j] = sin(2*M_PI*x_aux);
      pi[j] = sin(2*M_PI*x_aux);
      exact_phi[j] = (sin(2*M_PI*(x_aux-simul_time)) + sin(2*M_PI*(x_aux+simul_time)))/2 + (1/4*M_PI) *  (cos(2*M_PI * (x_aux - simul_time)) - cos(2*M_PI*(x_aux-simul_time)));

     //cout << j << "     " << setw(5) << phi[j] << "     " <<setw(5) << exact_phi[j] << setw(5) << "     " << abs(phi[j]-exact_phi[j]) << "     " << Deviation(phi[j],exact_phi[j]) << endl;
     l1_norm = l1_norm + abs(phi[j]-exact_phi[j]);
    }

    cout << "dx = " << dx_i << ", l1 = " << l1_norm << " (" << idx << "th iteration)" << endl;

    CONVERGENCE_FILE << dx_i << " "<< l1_norm << endl; 
  }
  CONVERGENCE_FILE.close();  
}


int main(){

  // define fields
  double* y = new double[2*N];
  double* phi = &y[0];
  double* pi = &y[N];
  double* exact_phi = new double[N];

  // impose initial values on phi and pi
  double x_aux; 
  /*
  for(int j=0; j<N; j++){
    x_aux = j*dx;
        // for a sinusoidal
    phi[j] = sin(2*M_PI*x_aux);
    pi[j] = sin(2*M_PI*x_aux);
	  exact_phi[j] = (sin(2*M_PI*(x_aux-simul_time))+sin(2*M_PI*(x_aux+simul_time)))/2 + (1/4*M_PI) *  (cos(2*M_PI * (x_aux - simul_time)) - cos(2*M_PI*(x_aux-simul_time)));
  */


        // for a Gaussian
  //double x0=0.5;
	//double sigma=5;
	//phi[j] = Gaussian(x_aux,sigma,x0);
	//pi[j] = 0;
	//exact_phi[j] = 0.5*(Gaussian(x_aux-simul_time,sigma,x0) + Gaussian(x_aux+simul_time,sigma,x0));
  //}

  for(int i=0; i<number_iterations; i++){
      RungeKutta(y, RHS);
      //Saving_Data_File(y,i);
  }

    Convergence_Plot(y, exact_phi, 1E-7, 1E-4, 1000);

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
  //cout << "j" << setw(5) << "phi" << setw(5) << "Xphi" << endl;
  //for(int j=0;j<N;j++){
  //cout << j << "     " << setw(5) << phi[j] << "     " <<setw(5) << exact_phi[j] << setw(5) << "     " <<abs(phi[j]-exact_phi[j]) << "     " << Deviation(phi[j],exact_phi[j])<< endl;
  //}
  // delete pointer arrays
  delete[] y;
    
  return 0;
  
}
