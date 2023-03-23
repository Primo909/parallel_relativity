#include "WaveEquationSolver1D.h"

double Gaussian(double x, double sigma, double x0){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}

WaveEquationSolver1D::WaveEquationSolver1D(){}

void WaveEquationSolver1D::Saving_Data_File(int N, double* state_vector, double* axis, int iteration){

  string string_aux=File_String+to_string(iteration+1)+".dat";
  
  fstream ITERATION_FILE;
  ITERATION_FILE.open(string_aux,ios::out);

  if (!ITERATION_FILE){                 
      cout<<"Error: File not created"<<endl;    
  }

  for(int j=0; j<N; j++) ITERATION_FILE<<axis[j]<<" "<<state_vector[j]<<endl;

  ITERATION_FILE.close();   

}

void WaveEquationSolver1D::SetInitialConditions(int N, double *y, double* InitialState){
    for(int j=0; j<2*N; j++) y[j]=InitialState[j];
}

void WaveEquationSolver1D::SecondDerivative(int N, double dx, double* field, double* second_derivative_vector){
    int number_ghosts=3;
    double* right_ghost_cells_field = new double[number_ghosts];
    double* left_ghost_cells_field = new double[number_ghosts];
    // assign values to ghost cells
    for(int j=0; j<number_ghosts; j++){
        right_ghost_cells_field[j] = field[j];
        left_ghost_cells_field[number_ghosts-1-j] = field[N-1-j];
    }
    // calculate the second derivative with five points (-2,-1,0,1,2)
    // formula from https://web.media.mit.edu/~crtaylor/calculator.html
    for(int j=0; j<N; j++){
        if(j==0) second_derivative_vector[j] = (-1*left_ghost_cells_field[number_ghosts-2]+16*left_ghost_cells_field[number_ghosts-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
        else if(j==1) second_derivative_vector[j] = (-1*left_ghost_cells_field[number_ghosts+1-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
        else if(j==N-1) second_derivative_vector[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*right_ghost_cells_field[0]-1*right_ghost_cells_field[1])/(12*1.0*dx*dx);
        else if(j==N-2) second_derivative_vector[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*right_ghost_cells_field[0])/(12*1.0*dx*dx);
        else second_derivative_vector[j] = (-1*field[j-2]+16*field[j-1]-30*field[j+0]+16*field[j+1]-1*field[j+2])/(12*1.0*dx*dx);
    }
    // delete the ghost cells
    delete[] left_ghost_cells_field, right_ghost_cells_field;

}

void WaveEquationSolver1D::RHS(int N, double dx, double* state_vector, double* rhs_vector){
    double* phi = &state_vector[0];
    double* pi = &state_vector[N];
    for(int j=0; j<N; j++) rhs_vector[j] = pi[j];
    SecondDerivative(N, dx, phi, &rhs_vector[N]);
    //for(int j=0; j<N; j++) cout<<phi[j]<<" "<<rhs_vector[j]<<endl;
    //for(int j=0; j<N; j++) cout<<pi[j]<<" "<<rhs_vector[N+j]/(2*M_PI*2*M_PI)<<endl;
}

void WaveEquationSolver1D::RuggeKutta(int N, double dx, double dt, double* state_vector){
    double* k1 = new double[2*N];
    double* k2 = new double[2*N];
    double* k3 = new double[2*N];
    double* k4 = new double[2*N];

    double* vector_temp = new double[2*N];

    RHS(N,dx,state_vector, k1); // Calc k1
    for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dt*k1[j]/2;
    RHS(N,dx,vector_temp, k2); // Calc k2
    for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dt*k2[j]/2; 
    RHS(N,dx,vector_temp, k3); // Calc k3
    for(int j=0; j<2*N; j++) vector_temp[j] = state_vector[j]+dt*k3[j];
    RHS(N,dx,vector_temp, k4);// Calc k4
    for(int j=0; j<2*N; j++) state_vector[j] = state_vector[j]+(k1[j] + 2*k2[j] + 2*k3[j] +k4[j])*dt/6; // Overwrite the state vector

    delete[] vector_temp,k1,k2,k3,k4;
}

void WaveEquationSolver1D::Solve(double xmin, double xmax, double dx, double* InitialState, double dt=1E-3, double simul_time=2){

    auto start = high_resolution_clock::now(); //Counting time
    cout<<Simulation_Start_String<<endl;

    int N=(xmax-xmin)/dx;
    int number_iterations = simul_time / dt;

    double* y= new double[2*N];
    double *phi= &y[0];
    double *pi= &y[N];

    double* axis= new double[N];
    for(int j=0; j<N; j++) axis[j]=xmin+j*dx;

    SetInitialConditions(N, y, InitialState);

    for(int i=0; i<number_iterations; i++){
        RuggeKutta(N,dx,dt,y);
        Saving_Data_File(N,y,axis,i);
    }

    auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout<<endl;
    cout<<Simulation_End_String<<":  "<<duration.count()/1E6<<" s"<<endl;

    delete[] axis,y;

}

void WaveEquationSolver1D::PointConvergenceTest(double (*func1)(double), double (*func2)(double), string filename){

    auto start = high_resolution_clock::now(); //Counting time
    cout<<Simulation_Start_String<<endl;
  
    fstream ITERATION_FILE;
    ITERATION_FILE.open(filename,ios::out);

    if (!ITERATION_FILE){
		cout<<"Error: File not created"<<endl;    
  }

	// T = simulation time
    double T=0.5;
    double dt=1E-4;
    int number_iterations= T/dt;

    double x_min=-1;
    double x_max=1;

    double dx_low=5E-3; //lowest resolution
    int N_low=(x_max-x_min)/dx_low;
    double dx_mid=dx_low/2; //middle resolution
    int N_mid=(x_max-x_min)/dx_mid;
    double dx_high=dx_mid/2; //highest resolution
    int N_high=(x_max-x_min)/dx_high;

    double* y_low= new double[2*N_low];
    double* y_mid= new double[2*N_mid];
    double* y_high = new double[2*N_high];


	double* phi_ratio = new double[N_low];
	double* phi_low_mid = new double[N_low];
	double* phi_mid_high = new double[N_low];

    double sigma=1;
    double x_aux;
    double x0=0;

    for(int j=0; j<N_low; j++){ 
        x_aux=x_min+j*dx_low;
        y_low[j]=func1(x_aux);
        y_low[j+N_low]=func2(x_aux);
    }
    for(int j=0; j<N_mid; j++){ 
        x_aux=x_min+j*dx_mid;
        y_mid[j]=func1(x_aux);
        y_mid[j+N_mid]=func2(x_aux);
    }    
    for(int j=0; j<N_high; j++){ 
        x_aux=x_min+j*dx_high;
        y_high[j]=func1(x_aux);
        y_high[j+N_high]=func2(x_aux);
    }

    for(int i=0; i<number_iterations; i++){
        RuggeKutta(N_low,dx_low,dt,y_low);
        RuggeKutta(N_mid,dx_mid,dt,y_mid);
        RuggeKutta(N_high,dx_high,dt,y_high);
    }

	for(int i=0; i<N_low; i++){
		phi_low_mid[i] = y_low[i] - y_mid[2*i];
		phi_mid_high[i] = y_mid[2*i] - y_high[4*i];
		phi_ratio[i] = phi_low_mid[i] / phi_mid_high[i];
		//cout << "     " << phi_low_mid[i] << "     " << phi_mid_high[i] << "    " << phi_ratio[i] << endl;
		ITERATION_FILE << "     " << phi_low_mid[i] << "     " << phi_mid_high[i] << "    " << phi_ratio[i] << endl;
	}

    auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout<<endl;
    cout<<Simulation_End_String<<":  "<<duration.count()/1E6<<" s"<<endl;

  ITERATION_FILE.close();   
}