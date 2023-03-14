#ifndef __GRID1D__
#define __GRID1D__

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class GRID1D{

    public:

    GRID1D(float, float, float);
    int NGrid();
    float dxGrid();

    private:

    float dx;
    float XMax;
    float XMin;
    int N;
    
}
;



#endif