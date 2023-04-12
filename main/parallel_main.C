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

double dx=1E-2; //Space discretization (uniform)
const double dt=1E-3; //Time discret.
const double x_min=-1, x_max=1; //Space interval 
const double simul_time = 2; //Simulation time
const int number_iterations = simul_time / dt; //Number of iterations
int N=(x_max - x_min) / dx; //Number of spatial cells
const int number_ghosts = 2; // Number of ghost cells
const double x0=(x_max+x_min)/2;
const double sigma=0.1;

double Sin(double x){
  return sin(2*M_PI*x);
}

double Zero(double x){
  return 0;
}

double GaussianFixed(double x){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}

void ParallelSecondDerivative(int n, double dx, double* field, double* second_derivative_field){

    int size,id;
    int my_first, my_last;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    //double* second_derivative_field=new double[n];

    double* guard_cells_send= new double[2*number_ghosts];
    double* guard_cells_recv= new double[2*number_ghosts];
    //Guard cell sending position [0 1][######][3 2]
    //Left guard cells
    guard_cells_send[0]=field[0];
    guard_cells_send[1]=field[1];
    //Right guard cells
    guard_cells_send[2]=field[n-2];
    guard_cells_send[3]=field[n-1];

    /*if(id==0){
        for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_send[j]<<endl;
    }*/

    int tag1=10,tag2=12;
    MPI_Request request1,request2;
    
    MPI_Status status1,status2;
    
    int left_id= (id-1+size)%size;
    int right_id= (id+1+size)%size;

    //cout<<left_id<<" "<<id<<" "<<right_id<<endl;

    ////AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    MPI_Isend(&guard_cells_send[0], number_ghosts, MPI_DOUBLE, left_id, tag1, MPI_COMM_WORLD, &request1);
    MPI_Irecv(&guard_cells_recv[0], number_ghosts, MPI_DOUBLE, right_id, tag1, MPI_COMM_WORLD, &request1);
    MPI_Isend(&guard_cells_send[2], number_ghosts, MPI_DOUBLE, right_id, tag2, MPI_COMM_WORLD, &request2);
    MPI_Irecv(&guard_cells_recv[2], number_ghosts, MPI_DOUBLE, left_id, tag2, MPI_COMM_WORLD, &request2);

    MPI_Wait(&request1, &status1);
    MPI_Wait(&request2, &status2);

    //for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_send[j]<<" "<<guard_cells_recv[j]<<endl;

    /*if(id==1){
        for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_recv[j]<<endl;
    }*/

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
    
    int size,id;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double* phi = &state_vector[0];
    double* pi = &state_vector[n];
    
    for(int j=0; j<n; j++) rhs_vector[j] = pi[j];
    //for(int j=0; j<n; j++) rhs_vector[j+n]= 1;
    ParallelSecondDerivative(n, dx, phi, &rhs_vector[n]);
    //for(int j=0; j<N; j++) cout<<phi[j]<<" "<<rhs_vector[j]<<endl;
    //for(int j=0; j<N; j++) cout<<pi[j]<<" "<<rhs_vector[N+j]/(2*M_PI*2*M_PI)<<endl;
    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized

}

void ParallelRungeKutta(int n, double dx, double dt, double* state_vector){
        
    int size,id;
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



int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);

    int size,id;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double startwtime = 0.0, endwtime;
    if(id==0){
        startwtime = MPI_Wtime();
        cout<<"Parallel "<< Simulation_Start_String<<endl;
    }

    int n=N/size;

    double* y= new double[2*n];
    double* phi= &y[0];
    double* pi= &y[n];

    double* axis= new double[n];
    for(int j=0; j<n; j++) {
        axis[j]=x_min+id*n*dx+j*dx;
    }

    for(int j=0; j<n; j++){
        y[j]=GaussianFixed(axis[j]); //initial condition phi
        y[n+j]=Zero(axis[j]); //initial condition pi
        //if(id==0) cout<<y[j]<<"  "<<y[j+N]<<endl;
    }

    MPI_Barrier(MPI_COMM_WORLD); //make sure they are syncronized
    
    double *PHI = new double[N];
    double *AXIS= new double[N];

    MPI_Gather(axis,n,MPI_DOUBLE,AXIS,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for(int i=0; i<number_iterations; i++){

        string string_aux=File_String+to_string(i+1)+".dat";

        fstream ITERATION_FILE;
        ITERATION_FILE.open(string_aux,ios::out);

        if (!ITERATION_FILE){                 
            cout<<"Error: File not created"<<endl;    
        }

        //if(id==1) cout<<y[10]<<endl;

        ParallelRungeKutta(n, dx, dt, y);
        MPI_Gather(y,n,MPI_DOUBLE,PHI,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(id==0){
            for(int j=0; j<N; j++) ITERATION_FILE<<AXIS[j]<<" "<<PHI[j]<<endl;
        }

        ITERATION_FILE.close();

    }



    if(id==0){
        endwtime = MPI_Wtime();
        cout<<Simulation_End_String<<endl;
        cout<<endwtime-startwtime<<endl;
    }

    /*double* second_derivative_phi=new double[n];

    ParallelSecondDerivative(n,dx,phi,second_derivative_phi);
    */

    /*double* guard_cells_send= new double[2*number_ghosts];
    double* guard_cells_recv= new double[2*number_ghosts];
    //Guard cell sending position [0 1][######][3 2]
    //Left guard cells
    guard_cells_send[0]=phi[0];
    guard_cells_send[1]=phi[1];
    //Right guard cells
    guard_cells_send[2]=phi[n-1];
    guard_cells_send[3]=phi[n-2];

    /*if(id==0){
        for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_send[j]<<endl;
    }

    int tag1=10,tag2=12;
    MPI_Request request1,request2;
    
    MPI_Status status1,status2;
    
    int left_id= (id-1+size)%size;
    int right_id= (id+1+size)%size;

    //cout<<left_id<<" "<<id<<" "<<right_id<<endl;

    MPI_Isend(&guard_cells_send[0], number_ghosts, MPI_DOUBLE, left_id, tag1, MPI_COMM_WORLD, &request1);
    MPI_Irecv(&guard_cells_recv[2], number_ghosts, MPI_DOUBLE, right_id, tag1, MPI_COMM_WORLD, &request1);
    MPI_Isend(&guard_cells_send[2], number_ghosts, MPI_DOUBLE, right_id, tag2, MPI_COMM_WORLD, &request2);
    MPI_Irecv(&guard_cells_recv[0], number_ghosts, MPI_DOUBLE, left_id, tag2, MPI_COMM_WORLD, &request2);

    MPI_Wait(&request1, &status1);
    MPI_Wait(&request2, &status2);

    //for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_send[j]<<" "<<guard_cells_recv[j]<<endl;

    /*if(id==1){
        for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_recv[j]<<endl;
    }

    for(int j=0; j<n; j++){
        if(j==0) second_derivative_phi[j] = (-1*guard_cells_recv[0]+16*guard_cells_recv[1]-30*phi[j+0]+16*phi[j+1]-1*phi[j+2])/(12*1.0*dx*dx);
        else if(j==1) second_derivative_phi[j] = (-1*guard_cells_recv[1]+16*phi[j-1]-30*phi[j+0]+16*phi[j+1]-1*phi[j+2])/(12*1.0*dx*dx);
        else if(j==n-1) second_derivative_phi[j] = (-1*phi[j-2]+16*phi[j-1]-30*phi[j+0]+16*guard_cells_recv[2]-1*guard_cells_recv[3])/(12*1.0*dx*dx);
        else if(j==n-2) second_derivative_phi[j] = (-1*phi[j-2]+16*phi[j-1]-30*phi[j+0]+16*phi[j+1]-1*guard_cells_recv[2])/(12*1.0*dx*dx);
        else second_derivative_phi[j] = (-1*phi[j-2]+16*phi[j-1]-30*phi[j+0]+16*phi[j+1]-1*phi[j+2])/(12*1.0*dx*dx);
    }
    */
    MPI_Finalize();

    return 0;
}