#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Var1D.h"
#include "Grid1D.h"

using namespace std;

float dx=1E-2;
float dt=1E-2;
float XMin=0;
float XMax=1;
int N=(XMax-XMin)/dx;


int main(){

  GRID1D GRID= GRID1D(XMax,XMin,dx);

  float* X= new float[N];
  float* dXdx= new float[N];
  float* ddXddx=new float[N];

  for(int j=0; j<N; j++) X[j]=sin(2*M_PI*j*dx);

  Var1D VarX= Var1D(GRID, X);
  VarX.SetGhostCells();
  VarX.Derivative(dXdx,2); 
  VarX.SecondDerivative(ddXddx,2);

  //for(int j=0; j<N; j++) cout<<X[j]<<" "<<2*M_PI*cos(2*M_PI*j*dx)<<" "<<dXdx[j]<<endl;
  for(int j=0; j<N; j++) cout<<X[j]<<" "<<-2*M_PI*2*M_PI*X[j]<<" "<<ddXddx[j]<<endl;

  delete[] X,dXdx,ddXddx;

  return 0;
  
}
