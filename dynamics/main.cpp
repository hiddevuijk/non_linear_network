// nr3 headers
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
#include "nr_headers/fourier.h"
#include "nr_headers/correl.h"
#include "nr_headers/stepper.h"
#include "nr_headers/stepperbs.h"
#include "nr_headers/odeint.h"
#include "nr_headers/stepperdopr853.h"
#include "nr_headers/stepperdopr5.h"
#include "nr_headers/stepperstoerm.h"

// other headers
#include "headers/derivatives.h"
#include "headers/write_matrix.h"
#include "headers/box_muller.h"
#include "headers/read_input.h"
#include "headers/functions.h"
#include "headers/other.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>

using namespace::std;


int main(int argc, char* argv[])
{

	int N;		// number of neurons
	double g;	// gain parameter
	int seed;	// seed for random generator

	int tf;		// integrate from 0 intil tf
	int tsave;	// number of point to save
	
	int function; // transfer function, defined int functions

	double mean_noise;
	double var_noise;

	// read valuese of the variables from input file
	read_input(N,g,tf,tsave,function,mean_noise,var_noise,seed, "input.txt");
	// check for commanline input
	if(argc ==2) {
		try {
			g = stod(argv[1]);
		} catch( exception& e) {
			cerr << "Unable to cast " << argv[1] << " to int" << endl;
			cerr << "\t error: " << e.what() << endl;
			return 1;
		}
	}

	double std_noise = sqrt(var_noise);
	double dt = tf/(double)tsave;
	double sqrt_dt = sqrt(dt);

	Ran r(seed);	// random number generator

	// standard deviation of connectivity matrix
	double std  = g/sqrt((double)N);

	// create connectivity matrix w
	vector<double> temp(N,.0);
	vector<vector<double> > w(N,temp);

	// state vector of the neurons
	VecDoub x(N,0.0);
	vector<double> tval(tsave,0.0);
	vector<vector<double> > xt(N,tval);

	// initialize w and x
	for(int i=0;i<N;++i){
		x[i] = r.doub(); // uniform distributed random number
		for( int j=0;j<N;++j) {
			w[i][j] =  std*bm_transform(r); // normal dist mean=0 std=std
		}
	}


	// integration settings
	double atol = 1.e-9;
	double rtol = atol;
	double h1 = 0.01;
	double hmin = 0.;

	// f is a pointer to the transfer function
	double (*f)(double);
	if(function == 1) f = th_linear;
	if(function == 2) f = tanh;
	if(function == 3) f = tanh01;

	NW nw(w,N,f);
	Output out;
	for(int ti=0;(ti+1)<tsave;++ti){
		Odeint<StepperDopr5<NW> > ode(x,ti*dt,(ti+1)*dt,atol,rtol,h1,hmin,out,nw);
		ode.integrate();
		for(int i=0;i<N;++i) {
			xt[i][ti] = x[i];
			x[i] += sqrt_dt*std_noise*bm_transform(r);
		}
		tval[ti+1]= ti*dt;
	}

	

	// save results
	string add_to_name = "";

	stringstream stream;
	if(argc ==2){
		stream << fixed << setprecision(2) << g;
		add_to_name = '_' + stream.str(); 
	}
	write_matrix(tval,tsave,"t" + add_to_name + ".csv");
	write_matrix(xt,N,tsave,"x" + add_to_name + ".csv");
	
	return 0;
}


