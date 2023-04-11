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

void SecondDerivative(int n, double dx, double* field, double* second_derivative_field){

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
    guard_cells_send[2]=field[n-1];
    guard_cells_send[3]=field[n-2];

    /*if(id==0){
        for(int j=0; j<2*number_ghosts; j++) cout<<id<<" "<<guard_cells_send[j]<<endl;
    }*/

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
    }*/

    for(int j=0; j<n; j++){
        if(j==0) second_derivative_field[j] = (-1*guard_cells_recv[0]+16*guard_cells_recv[1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
        else if(j==1) second_derivative_field[j] = (-1*guard_cells_recv[1]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
        else if(j==n-1) second_derivative_field[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*guard_cells_recv[2]-1*guard_cells_recv[3])/(12*1.0*dx*dx);
        else if(j==n-2) second_derivative_field[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*guard_cells_recv[2])/(12*1.0*dx*dx);
        else second_derivative_field[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
    }
    
    // delete the ghost cells
    delete[] guard_cells_recv, guard_cells_send;

}

int main(int argc, char* argv[]){

//This part keep sequential for now, maybe change later
//
    MPI_Init(&argc, &argv);

    int size,id;
    int my_first, my_last;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    int n=N/size;

    my_first = id*n;
    my_last = (id+1)*n;

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

    double* second_derivative_phi=new double[n];

    SecondDerivative(n,dx,phi,second_derivative_phi);

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