#include "Grid1D.h"

GRID1D::GRID1D(float XMAX, float XMIN, float DX){

    XMax=XMAX;
    XMin=XMIN;
    dx=DX;
    N=(XMax-XMin)/dx;

}

int GRID1D::NGrid(){return N;}

float GRID1D::dxGrid(){return dx;}