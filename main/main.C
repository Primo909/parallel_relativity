#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
using namespace std;

double dx=1E-2;
double dt=1E-3;
double x_min=0;
double x_max=1;
double simul_time = 0.1;
int number_iterations = simul_time / dt;
int N=(x_max - x_min) / dx;
int number_ghosts = 3; // Number of ghost cells


int main(){

  // define fields
  double* phi = new double[N];
  double* pi = new double[N];

  // create ghost cells
  double* left_ghost_cells_phi = new double[number_ghosts];
  double* left_ghost_cells_pi = new double[number_ghosts];
  double* right_ghost_cells_phi = new double[number_ghosts];
  double* right_ghost_cells_pi = new double[number_ghosts];

  // impose initial values on phi and pi
  double x_aux;
  for(int j=0; j<N; j++){
    x_aux = j*dx;
    phi[j] = sin(2*M_PI*x_aux);
    pi[j] = 2*M_PI*cos(2*M_PI*x_aux);
  }

  //for(int j=0;j<N;j++){
  //cout << j << "," << phi[j] << "," << pi[j] << endl;
  //}

  // set values of ghost cells
  for(int j=0; j<number_ghosts; j++){
    left_ghost_cells_phi[number_ghosts-1-j] = phi[N-1-j];
    left_ghost_cells_pi[number_ghosts-1-j] = phi[N-1-j];
    right_ghost_cells_phi[j] = phi[j];
    right_ghost_cells_pi[j] = pi[j];
  }

  // define state vector
  vector<double*>y;
  y.push_back(pi);
  y.push_back(phi);

  
  // Runge kutta
  double* k1_phi = new double[N];
  double* k1_pi = new double[N];
  vector<double*>k1;
  k1.push_back(k1_pi);
  k1.push_back(k1_phi);

  double* k2_phi = new double[N];
  double* k2_pi = new double[N];
  vector<double*>k2;
  k2.push_back(k2_pi);
  k2.push_back(k2_phi);

  double* k3_phi = new double[N];
  double* k3_pi = new double[N];
  vector<double*>k3;
  k3.push_back(k3_pi);
  k3.push_back(k3_phi);

  double* k4_phi = new double[N];
  double* k4_pi = new double[N];
  vector<double*>k4;
  k4.push_back(k4_pi);
  k4.push_back(k4_phi);

  for(int i=0; i<number_iterations; i++){
    for(int j=0; j<N; j++){
      k1[0][j] = phi[j];
      k1[1][j] = pi[j];

      k2[0][j] = phi[j] + dt * k1[1][j]/2;
      k2[1][j] = pi[j] + dt * k1[0][j]/2;

      k3[0][j] = phi[j] + dt * k2[1][j]/2;
      k3[1][j] = pi[j] + dt * k2[0][j]/2;

      k4[0][j] = phi[j] + dt * k3[1][j];
      k4[1][j] = pi[j] + dt * k3[0][j];

    }
    for(int j=0; j<N; j++){
      y[0][j] = y[0][j] + (k1[0][j] + 2*k2[0][j] + 2*k3[0][j] + k4[0][j])*dt/6;
      y[1][j] = y[1][j] + (k1[1][j] + 2*k2[1][j] + 2*k3[1][j] + k4[1][j])*dt/6;
    }
  }
  // output something
  cout << "j" << "," << "phi" << "," << "pi" << endl;
  for(int j=0;j<N;j++){
  cout << j << "," << phi[j] << "," << pi[j] << endl;
  }
  // delete pointer arrays
  delete[] phi, pi;
  delete[] left_ghost_cells_phi, left_ghost_cells_pi, right_ghost_cells_phi, right_ghost_cells_pi;
  
  return 0;
  
}
