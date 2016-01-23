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

#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>

using namespace::std;

int main()
{

	int N;		// number of neurons
	double g;	// gain parameter
	int seed;	// seed for random generator

	int tf;		// integrate from 0 intil tf
	int tsave;	// number of point to save
	
	// read valuese of the variables from input file
	read_input(N,g,tf,tsave,seed, "input.txt");
	Ran r(seed);	// random number generator
	double std  = g/sqrt((double)N);	//std of connectivity matrix

	// create connectivity matrix w
	vector<double> temp(N,.0);
	vector<vector<double> > w(N,temp);

	// state vector of the neurons
	VecDoub x(N,0.0);

	// initialize w and x
	for(int i=0;i<N;++i){
		x[i] = r.doub(); // uniform distributed random number
		for( int j=0;j<N;++j) {
			w[i][j] = std*bm_transform(r); // normal dist mean=0 std=std
		}
	}


	// integration settings
	double atol = 1.e-9;
	double rtol = atol;
	double h1 = 0.01;
	double hmin = 0.;

	
	NW nw(w,N);
	Output out(tsave);
	Odeint<StepperDopr5<NW> > ode(x,0,tf,atol,rtol,h1,hmin,out,nw);
	ode.integrate();

	write_matrix(out.xsave,out.count,"t.csv");
	write_matrix(out.ysave,N,out.count,"x.csv");

	return 0;
}


