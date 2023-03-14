#include "Var1D.h"


Var1D::Var1D(GRID1D GRID, float* v){
    N=GRID.NGrid();
    dx=GRID.dxGrid();
    V= new float[N];
    dVdx= new float[N];
    ddVddx= new float[N];
    for(int j=0; j<N; j++) V[j]=v[j];
}

void Var1D::SetGhostCells(){
    RightGhostCells=new float[NGhostCells];
    LeftGhostCells=new float[NGhostCells];
    for(int j=0; j<NGhostCells; j++){
        LeftGhostCells[j]=V[N-NGhostCells+j];
        RightGhostCells[j]=V[j];
    }
}

void Var1D::Derivative(float* v, int O){
    if(O==1){
        for(int j=0; j<N; j++){
            if(j!=0 && j!=N-1) dVdx[j]=(V[j+1]-V[j-1])/(2*dx);
            else if(j==0) dVdx[j]= (V[j+1]-LeftGhostCells[NGhostCells-1])/(2*dx);
            else if(j==N-1) dVdx[j]= (RightGhostCells[0]-V[j-1])/(2*dx);
            v[j]=dVdx[j];
        }
    }else if(O==2){
        for(int j=0; j<N; j++){
            if(j!=0 && j!=1 && j!=N-1 && j!=N-2) dVdx[j]=(-V[j+2]+8*V[j+1]-8*V[j-1]+V[j-2])/(12*dx);
            else if(j==0) dVdx[j]=(-V[j+2]+8*V[j+1]-8*RightGhostCells[NGhostCells-1]+RightGhostCells[NGhostCells-2])/(12*dx);
            else if(j==1) dVdx[j]=(-V[j+2]+8*V[j+1]-8*V[j-1]+RightGhostCells[NGhostCells-1])/(12*dx);
            else if(j==N-2) dVdx[j]=(-RightGhostCells[0]+8*V[j+1]-8*V[j-1]+V[j-2])/(12*dx);
            else if(j==N-1) dVdx[j]=(-RightGhostCells[1]+8*RightGhostCells[0]-8*V[j-1]+V[j-2])/(12*dx);
            v[j]=dVdx[j];
        }
    }
}

void Var1D::SecondDerivative(float*v, int O){
    if(O==1){
        for(int j=0; j<N; j++){
            if(j!=0 && j!=N-1) ddVddx[j]=(V[j+1]-2*V[j]+V[j-1])/(dx*dx);
            else if(j==0) ddVddx[j]= (V[j+1]-2*V[j]+LeftGhostCells[NGhostCells-1])/(dx*dx);
            else if(j==N-1) ddVddx[j]= (RightGhostCells[0]-2*V[j]+V[j-1])/(dx*dx);
            v[j]=ddVddx[j];
        }
    }else if(O==2){
        for(int j=0; j<N; j++){
            if(j!=0 && j!=1 && j!=N-1 && j!=N-2) ddVddx[j]=(-V[j+2]+16*V[j+1]-30*V[j]+16*V[j-1]-V[j-2])/(12*dx*dx);
            else if(j==0) ddVddx[j]=(-V[j+2]+16*V[j+1]-30*V[j]+16*RightGhostCells[NGhostCells-1]-RightGhostCells[NGhostCells-2])/(12*dx*dx);
            else if(j==1) ddVddx[j]=(-V[j+2]+16*V[j+1]-30*V[j]+16*V[j-1]-RightGhostCells[NGhostCells-1])/(12*dx*dx);
            else if(j==N-2) ddVddx[j]=(-RightGhostCells[0]+16*V[j+1]-30*V[j]+16*V[j-1]-V[j-2])/(12*dx*dx);
            else if(j==N-1) ddVddx[j]=(-RightGhostCells[1]+16*RightGhostCells[0]-30*V[j]+16*V[j-1]-V[j-2])/(12*dx*dx);
            v[j]=ddVddx[j];
        }
    }
}