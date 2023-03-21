#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
#include <string>

using namespace std;

const double dx=1E-3; //Space discretization (uniform)
const double dt=1E-3; //Time discret.
const double x_min=0, x_max=1; //Space interval 
const double simul_time = 0.1; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
const int N=(x_max - x_min) / dx; //Number of spatial cells
const int number_ghosts = 3; // Number of ghost cells


string File_String= "numerical_phi_ite";

void Saving_Data_File(double* state_vector, int iteration){

    string iteration_s = to_string(iteration);
    string string_aux = File_String + iteration_s;
    
    fstream ITERATION_FILE;
    ITERATION_FILE.open(string_aux,ios::out);

    for(int j=0; j<N; j++) ITERATION_FILE << j*dx << " " << state_vector[j] << endl;

    ITERATION_FILE.close();   

}

int main(){

    return 0;
}
