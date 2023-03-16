#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip> 
using namespace std;

int N=3;

int main(){

  // define fields
  double* phi = new double[N];
  double* pi = new double[N];
  double* y = new double[2*N];
 
  for(int j=0; j<N; j++){
    phi[j] = 1;
    pi[j] = 2;
  }
  for(int j=0; j<N; j++){
    y[j] = 1;
    y[N+j] = 2;
  }

  double* y1 = &y[0] ;
  double* y2 = &y[N] ;

  cout << " " << "," << "phi" << "," << "pi" << endl;
  for(int j=0;j<N;j++){
  cout << j << "," << phi[j] << "," << pi[j] << endl;
  }

  cout << endl;

  cout << " " << "," << "y1" << "," << "y2" << endl;
  for(int j=0;j<N;j++){
  cout << j << "," << y1[j] << "," << y2[j] << endl;
  }
  
  
  cout << endl;
  cout << " ," << y << endl;
  for(int j=0;j<2*N;j++){
  cout << j << "," << y[j] << endl;
  }
  
  y[3] = 9;
  y1[1] = 66;
  cout << endl;
  cout << " ,y" << endl;
  for(int j=0;j<2*N;j++){
  cout << j << "," << y[j] << endl;
  }
  
  delete[] phi, pi;

  
  return 0;
  
}
