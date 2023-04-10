#include "PDESolver1D.h"

/*double Gaussian(double x, double sigma, double x0){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-0.5 * (x-x0)*(x-x0)/sigma/sigma);
}*/

PDESolver1D::PDESolver1D(double DX, double XMIN, double XMAX, double(*func3)(double ), double(*func4)(double),void(*func1)(int, double*,double*), void(*func2)(int, double*,double*)){

    dx=DX;    
    x_min=XMIN;
    x_max=XMAX;
    N=(x_min-x_max)/dx;
    RHS1=func1;
    RHS2=func2;
    Initial_Condition_Phi=func3;
    Initial_Condition_Pi=func4;
    Axis= new double[N];
    y= new double[2*N];
    //for(int j=0; j<N; j++) Axis[j]= x_min+j*dx;
}

PDESolver1D::~PDESolver1D(){
    //delete[] Axis;
    //delete[] y;
}

/*void PDESolver1D::Info_Conv_Test(double dx, string file_name){
    cout<<endl;
    cout<<Convergence_Test_String<<endl;
    cout<<"dt= 1E-4 T=1.2"<<endl;
    cout<<"Results printed on "<<file_name<<endl;
}*/

void PDESolver1D::Saving_Data_File(int iteration){

  string string_aux=File_String+to_string(iteration+1)+".dat";
  
  fstream ITERATION_FILE;
  ITERATION_FILE.open(string_aux,ios::out);

  if (!ITERATION_FILE){                 
      cout<<"Error: File not created"<<endl;    
  }

  for(int j=0; j<N; j++) ITERATION_FILE<<Axis[j]<<" "<<y[j]<<endl;

  ITERATION_FILE.close();   

}

void PDESolver1D::SetInitialConditions(){
    for(int i=0; i<N; i++){
        y[i]=Initial_Condition_Phi(Axis[i]);
        y[N+i]=Initial_Condition_Pi(Axis[i]);
    }    
}

void PDESolver1D::RHS(double* state_vector, double* rhs_vector){
    double* phi = &state_vector[0];
    double* pi = &state_vector[N];
    
    RHS1(N, phi, &rhs_vector[0]);
    RHS2(N, pi, &rhs_vector[N]);
}

void PDESolver1D::RuggeKutta(double dt){
    double* k1 = new double[2*N];
    double* k2 = new double[2*N];
    double* k3 = new double[2*N];
    double* k4 = new double[2*N];

    double* vector_temp = new double[2*N];

    RHS(y, k1); // Calc k1
    for(int j=0; j<2*N; j++) vector_temp[j] = y[j]+dt*k1[j]/2;
    RHS(vector_temp, k2); // Calc k2
    for(int j=0; j<2*N; j++) vector_temp[j] = y[j]+dt*k2[j]/2; 
    RHS(vector_temp, k3); // Calc k3
    for(int j=0; j<2*N; j++) vector_temp[j] = y[j]+dt*k3[j];
    RHS(vector_temp, k4);// Calc k4
    for(int j=0; j<2*N; j++) y[j] = y[j]+(k1[j] + 2*k2[j] + 2*k3[j] +k4[j])*dt/6; // Overwrite the state vector

    delete[] vector_temp,k1,k2,k3,k4;
}

void PDESolver1D::Solve(double dt=1E-3, double simul_time=2, bool file=false){

    auto start = high_resolution_clock::now(); //Counting time
    cout<<Simulation_Start_String<<endl;

    int number_iterations = simul_time / dt;

    //SetInitialConditions();

    /*for(int i=0; i<number_iterations; i++){
        RuggeKutta(dt);
        if(file==1) Saving_Data_File(i);
    }*/

    auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout<<endl;
    cout<<Simulation_End_String<<":  "<<duration.count()/1E6<<" s"<<endl;

}

/*void PDESolver1D::PointConvergenceTest(double (*func1)(double), double (*func2)(double), string filename){

    auto start = high_resolution_clock::now(); //Counting time
    cout<<Simulation_Start_String<<endl;
  
    fstream POINT_CONV_FILE;
    POINT_CONV_FILE.open(filename,ios::out);

    if (!POINT_CONV_FILE){
		cout<<"Error: File not created"<<endl;    
    }

	// T = simulation time
    double T=1.2;
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

	double phi_ratio; 
	double phi_low_mid;
	double phi_mid_high;

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
		phi_low_mid = y_low[i] - y_mid[2*i];
		phi_mid_high = y_mid[2*i] - y_high[4*i];
		phi_ratio = phi_low_mid / phi_mid_high;
		POINT_CONV_FILE << "     " << phi_low_mid << "     " << phi_mid_high << "    " << phi_ratio << endl;
	}

    POINT_CONV_FILE.close();

    Info_Conv_Test(dx_low,filename);

    auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
    cout<<Simulation_End_String<<duration.count()/1E6<<" s"<<endl;

    delete[] y_low, y_mid, y_high;

}

void PDESolver1D::NormConvergenceTest(double (*func1)(double), double (*func2)(double), string filename){

    auto start = high_resolution_clock::now(); //Counting time
    fstream NORM_CONV_FILE;
    NORM_CONV_FILE.open(filename,ios::out);
    if (!NORM_CONV_FILE){
		cout<<"Error: File not created"<<endl;    
    }

	// T = simulation time
    double T=1.2;
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
    
    double phi_low_mid, phi_mid_high, sum_low_mid=0, sum_mid_high=0;
    double p;

    for(int j=0; j<number_iterations; j++){
        RuggeKutta(N_low,dx_low,dt,y_low);
        RuggeKutta(N_mid,dx_mid,dt,y_mid);
        RuggeKutta(N_high,dx_high,dt,y_high);

        sum_low_mid=0, sum_mid_high=0;

        for(int i=0; i<N_low; i++){
            phi_low_mid = (y_low[i] - y_mid[2*i])*(y_low[i] - y_mid[2*i]);
            phi_mid_high = (y_mid[2*i] - y_high[4*i])*(y_mid[2*i] - y_high[4*i]);
            sum_low_mid+=phi_low_mid;
            sum_mid_high+=phi_mid_high;
        }
        p=log2(sqrt(sum_low_mid)/sqrt(sum_mid_high));

        NORM_CONV_FILE<<j*dt<<" "<<p<<endl;

    }

    NORM_CONV_FILE.close();

    Info_Conv_Test(dx_low,filename);

    auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
    cout<<Simulation_End_String<<duration.count()/1E6<<" s"<<endl;

    delete[] y_low, y_mid, y_high;


}*/
