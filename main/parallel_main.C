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

using namespace std;
using namespace std::chrono;

int N = 1000; //Space discretization (uniform)
const double dt = 1E-4; //Time discret.
const double x_min = -1, x_max = 1; //Space interval 
const double simul_time = 2; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
double dx = 0;
const int number_ghosts = 2; // Number of ghost cells
const double x0 = (x_max+x_min)/2;
const double sigma = 0.1;


double Sin(double x){
  return sin(2*M_PI*x);
}


double Zero(double x){
  return 0;
}


double GaussianFixed(double x){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}


void Saving_Data_File(int N, double* state_vector, double* Axis, int iteration){

    string string_aux=File_String+to_string(iteration+1)+".dat";

    fstream ITERATION_FILE;
    ITERATION_FILE.open(string_aux,ios::out);

    if (!ITERATION_FILE){                 
        cout<<"Error: File not created"<<endl;    
    }

    for(int j=0; j<N; j++) ITERATION_FILE << Axis[j] << " "<<state_vector[j] << endl;

    ITERATION_FILE.close();   
}


void ParallelSecondDerivative(int n, double dx, double* field, double* second_derivative_field){

    int size,id;
    int my_first, my_last;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double* guard_cells_send = new double[2*number_ghosts];
    double* guard_cells_recv = new double[2*number_ghosts];
    
    //Guard cell sending position [0 1][######][3 2]
    //Left guard cells
    guard_cells_send[0] = field[0];
    guard_cells_send[1] = field[1];
    //Right guard cells
    guard_cells_send[2] = field[n-2];
    guard_cells_send[3] = field[n-1];

    int tag1 =10,tag2 = 12;
    MPI_Request request1,request2;
    MPI_Status status1,status2;
    
    int left_id = (id - 1 + size)%size;
    int right_id = (id + 1 + size)%size;

    if(size>1){

        MPI_Isend(&guard_cells_send[0], number_ghosts, MPI_DOUBLE, left_id, tag1, MPI_COMM_WORLD, &request1);
        MPI_Irecv(&guard_cells_recv[0], number_ghosts, MPI_DOUBLE, right_id, tag1, MPI_COMM_WORLD, &request1);
        MPI_Isend(&guard_cells_send[2], number_ghosts, MPI_DOUBLE, right_id, tag2, MPI_COMM_WORLD, &request2);
        MPI_Irecv(&guard_cells_recv[2], number_ghosts, MPI_DOUBLE, left_id, tag2, MPI_COMM_WORLD, &request2);

        MPI_Wait(&request1, &status1);
        MPI_Wait(&request2, &status2);
    
    }
    
    for(int j=0; j<n; j++){
        if(j==0) second_derivative_field[j] = (-1*guard_cells_recv[2]+16*guard_cells_recv[3]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
        else if(j==1) second_derivative_field[j] = (-1*guard_cells_recv[3]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
        else if(j==n-1) second_derivative_field[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*guard_cells_recv[0]-1*guard_cells_recv[1])/(12*1.0*dx*dx);
        else if(j==n-2) second_derivative_field[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*guard_cells_recv[0])/(12*1.0*dx*dx);
        else second_derivative_field[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
    }
    
    // delete the ghost cells
    delete[] guard_cells_recv, guard_cells_send;
}


void ParallelRHS(int n, double dx, double* state_vector, double* rhs_vector){
    
    int size, id;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double* phi = &state_vector[0];
    double* pi = &state_vector[n];
    
    for(int j=0; j<n; j++) rhs_vector[j] = pi[j];
    ParallelSecondDerivative(n, dx, phi, &rhs_vector[n]);
    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized

}

void ParallelRungeKutta(int n, double dx, double dt, double* state_vector){
        
    int size, id;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double* k1 = new double[2*n];
    double* k2 = new double[2*n];
    double* k3 = new double[2*n];
    double* k4 = new double[2*n];

    double* vector_temp = new double[2*n];

    ParallelRHS(n,dx,state_vector, k1); // Calc k1
    for(int j=0; j<2*n; j++) vector_temp[j] = state_vector[j]+dt*k1[j]/2;
    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized
    ParallelRHS(n,dx,vector_temp, k2); // Calc k2
    for(int j=0; j<2*n; j++) vector_temp[j] = state_vector[j]+dt*k2[j]/2; 
    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized
    ParallelRHS(n,dx,vector_temp, k3); // Calc k3
    for(int j=0; j<2*n; j++) vector_temp[j] = state_vector[j]+dt*k3[j];
    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized
    ParallelRHS(n,dx,vector_temp, k4);// Calc k4
    for(int j=0; j<2*n; j++) state_vector[j] = state_vector[j]+(k1[j] + 2*k2[j] + 2*k3[j] +k4[j])*dt/6; // Overwrite the state vector
    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized

    delete[] vector_temp,k1,k2,k3,k4;
}


void PointConvergenceTest(string filename){

    fstream POINT_CONV_FILE;
    POINT_CONV_FILE.open(filename,ios::out);

    if (!POINT_CONV_FILE){
		cout<<"Error: File not created"<<endl;    
    }

    int size,id;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double T = 1.2;  // simulation time
    double dt = 1E-4;
    int number_iterations = T/dt;


    double dx_low = 5E-3;  // lowest resolution
    int N_low = (x_max - x_min)/dx_low;
    N_low=N_low-N_low%size;
    int n_low = N_low/size;
    double dx_mid = dx_low/2;  // middle resolution
    int N_mid = (x_max - x_min)/dx_mid;
    N_mid=N_mid-N_mid%size;
    int n_mid = N_mid/size;
    double dx_high = dx_mid/2;  // highest resolution
    int N_high =(x_max - x_min)/dx_high;
    N_high=N_high-N_high%size;
    int n_high = N_high/size;

    double* y_low = new double[2*n_low];
    double* y_mid = new double[2*n_mid];
    double* y_high = new double[2*n_high];
    double* axis = new double[n_low];

	double phi_ratio; 
	double phi_low_mid;
	double phi_mid_high;

    double x_aux;

    for(int j=0; j<n_low; j++){ 
        x_aux = x_min + id*n_low*dx_low + j*dx_low;
	axis[j]=x_aux;
        y_low[j] = GaussianFixed(x_aux);  // initial condition phi
        y_low[n_low + j] = Zero(x_aux);  // initial condition pi
    }
    for(int j=0; j<n_mid; j++){ 
        x_aux = x_min + id*n_mid*dx_mid + j*dx_mid;
        y_mid[j] = GaussianFixed(x_aux);  // initial condition phi
        y_mid[n_mid + j] = Zero(x_aux);  // initial condition pi
    } 
    for(int j=0; j<n_high; j++){ 
        x_aux = x_min + id*n_high*dx_high + j*dx_high;
        y_high[j] = GaussianFixed(x_aux);  // initial condition phi
        y_high[n_high + j] = Zero(x_aux);  // initial condition pi
    }

    MPI_Barrier(MPI_COMM_WORLD);  // make sure processes are syncronized

    for(int i=0; i<number_iterations; i++){
        ParallelRungeKutta(n_low, dx_low, dt, y_low);
        ParallelRungeKutta(n_mid, dx_mid, dt, y_mid);
        ParallelRungeKutta(n_high, dx_high, dt, y_high);
    }

    double *PHI_LOW = new double[N_low];
    double *PHI_MID = new double[N_mid];
    double *PHI_HIGH = new double[N_high];
    double *AXIS = new double[N_low];
    
    MPI_Gather(y_low, n_low, MPI_DOUBLE, PHI_LOW, n_low, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(y_mid, n_mid, MPI_DOUBLE, PHI_MID, n_mid, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(y_high, n_high, MPI_DOUBLE, PHI_HIGH, n_high, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(axis, n_low, MPI_DOUBLE, AXIS, n_low, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(id==0){
        for(int i=0; i<N_low; i++){
            phi_low_mid = PHI_LOW[i] - PHI_MID[2*i];
            phi_mid_high = PHI_MID[2*i] - PHI_HIGH[4*i];
            phi_ratio = phi_low_mid / phi_mid_high;
            POINT_CONV_FILE << AXIS[i] << "     " << phi_low_mid << "     " << phi_mid_high << "    " << phi_ratio << endl;
        }
    }

    POINT_CONV_FILE.close();

    delete[] y_low, y_mid, y_high, PHI_HIGH, PHI_LOW, PHI_MID, axis, AXIS;
}


void NormConvergenceTest(string filename){

    fstream NORM_CONV_FILE;
    NORM_CONV_FILE.open(filename,ios::out);

    if (!NORM_CONV_FILE){
		cout<<"Error: File not created"<<endl;    
    }

    int size, id;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double T = 1.2;  // simulation time
    double dt = 1E-4;
    int number_iterations = T/dt;

    double dx_low = 5E-3;  // lowest resolution
    int N_low = (x_max - x_min)/dx_low;
    int n_low = N_low/size;
    if(N_low%size!=0 && id==0) n_low = n_low+N_low%size;
    double dx_mid = dx_low/2;  // middle resolution
    int N_mid = (x_max - x_min)/dx_mid;
    int n_mid = N_mid/size;
    if(N_mid%size!=0 && id==0) n_mid = n_mid+N_mid%size;
    double dx_high = dx_mid/2;  // highest resolution
    int N_high = (x_max - x_min)/dx_high;
    int n_high = N_high/size;
    if(N_high%size!=0 && id==0) n_high = n_high+N_high%size;

    double* y_low = new double[2*n_low];
    double* y_mid = new double[2*n_mid];
    double* y_high = new double[2*n_high];

	double phi_ratio; 
	double phi_low_mid;
	double phi_mid_high;

    double sum_low_mid = 0, sum_mid_high = 0;
    double p = 0;

    double x_aux;

    for(int j=0; j<n_low; j++){ 
        x_aux = x_min + id*n_low*dx_low + j*dx_low;
        y_low[j] = GaussianFixed(x_aux);  // initial condition phi
        y_low[n_low + j] = Zero(x_aux);  // initial condition pi
    }
    for(int j=0; j<n_mid; j++){ 
        x_aux = x_min + id*n_mid*dx_mid + j*dx_mid;
        y_mid[j] = GaussianFixed(x_aux);  // initial condition phi
        y_mid[n_mid + j] = Zero(x_aux);  // initial condition pi
    } 
    for(int j=0; j<n_high; j++){ 
        x_aux = x_min + id*n_high*dx_high + j*dx_high;
        y_high[j] = GaussianFixed(x_aux);  // initial condition phi
        y_high[n_high + j] = Zero(x_aux);  // initial condition pi
    }

    MPI_Barrier(MPI_COMM_WORLD);  // make sure processes are syncronized

    double *PHI_LOW = new double[N_low];
    double *PHI_MID = new double[N_mid];
    double *PHI_HIGH = new double[N_high];

    for(int i=0; i<number_iterations; i++){
        ParallelRungeKutta(n_low, dx_low ,dt, y_low);
        ParallelRungeKutta(n_mid, dx_mid, dt, y_mid);
        ParallelRungeKutta(n_high, dx_high, dt, y_high);

        MPI_Gather(y_low, n_low, MPI_DOUBLE, PHI_LOW, n_low, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(y_mid, n_mid, MPI_DOUBLE, PHI_MID, n_mid, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(y_high, n_high, MPI_DOUBLE, PHI_HIGH, n_high, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(id==0){

            sum_low_mid = 0, sum_mid_high = 0;
            
            for(int i=0; i<N_low; i++){
                phi_low_mid = (PHI_LOW[i] - PHI_MID[2*i])*(PHI_LOW[i] - PHI_MID[2*i]);
                phi_mid_high = (PHI_MID[2*i] - PHI_HIGH[4*i])*(PHI_MID[2*i] - PHI_HIGH[4*i]);
                sum_low_mid += phi_low_mid;
                sum_mid_high += phi_mid_high;
            }

            p = log2(sqrt(sum_low_mid)/sqrt(sum_mid_high));
            NORM_CONV_FILE << i*dt << "  " << p << endl;
        }
    }

    NORM_CONV_FILE.close();

    delete[] y_low, y_mid, y_high, PHI_HIGH, PHI_LOW, PHI_MID;
}

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);

    int size,id;

    int N = atoi(argv[1]);
    bool saving = atoi(argv[2]);
    bool point_conv_test = atoi(argv[3]);
    double dx=(x_max-x_min)/N;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(id==0) cout << "Computing on " << size << " cores." << endl;
    if(id==0) cout << "dx = " << dx << endl;
    double startwtime = 0.0, endwtime;

    if(id==0){
        startwtime = MPI_Wtime();
        cout << "Parallel " << Simulation_Start_String << endl;
    }

    N = N - N%size;
    int n = N/size;
    cout << id << " " << n << endl;


    double* y = new double[2*n];
    double* phi = &y[0];
    double* pi = &y[n];

    double* axis = new double[n];
    for(int j=0; j<n; j++) {
        axis[j] = x_min + id*n*dx + j*dx;
    }

    for(int j=0; j<n; j++){
        y[j] = GaussianFixed(axis[j]);  // initial condition phi
        y[n+j] = Zero(axis[j]);  // initial condition pi
    }

    MPI_Barrier(MPI_COMM_WORLD);  //make sure processes are syncronized
    
    double *PHI = new double[N];
    double *AXIS = new double[N];

    MPI_Gather(axis, n, MPI_DOUBLE, AXIS, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int i=0; i<number_iterations; i++){

        ParallelRungeKutta(n, dx, dt, y);
        MPI_Gather(y, n, MPI_DOUBLE, PHI, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	    if(id==0 && i%100==0 && saving==true){
            //MPI_Gather(y, n, MPI_DOUBLE, PHI, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            Saving_Data_File(N,PHI,AXIS,i);
        }
    }

    if(id==0){
        endwtime = MPI_Wtime();
        cout << Simulation_End_String << endl;
        cout << "control," << atoi(argv[1]) << "," << size << "," << endwtime - startwtime<<endl;
    }

    if(point_conv_test==true){PointConvergenceTest("Data/pointConvTest.dat");}
    MPI_Finalize();

    return 0;
}
//  :) If you're reading this, you're awesome :) 
