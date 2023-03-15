#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
using namespace std;

double dx=1E-1;
double dt=1E-1;
double x_min=0;
double x_max=1;
double simul_time = 0.1;
int number_iterations = simul_time / dt;
int N=(x_max - x_min) / dx;
int number_ghosts = 3; // Number of ghost cells

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
  
  cout << "j   " << "phi   " << "pi" << endl;
  for(int j=0; j<N; j++){
  cout << j << "   " << phi[j] << "   " <<  pi[j] << endl;
  }
  
  cout << endl;
  
  y[4] = 66;
  phi[1] = 33;
  pi[2] = 55;
  cout << "j   " << "phi   " << "pi" << endl;
  for(int j=0; j<N; j++){
  cout << j << "   " << phi[j] << "   " <<  pi[j] << endl;
  }
  
  delete[] y, phi, pi;

  
  return 0;
  
}
