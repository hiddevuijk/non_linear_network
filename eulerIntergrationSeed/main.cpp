#include "Eigen/Core"

// numerical recipies headers (used for analysis)
#include "nr_headers/nr3.h"
#include "nr_headers/fourier.h"

// other headers
//#include "headers/derivatives.h"
#include "headers/write_matrix.h"
#include "headers/read_input.h"
#include "headers/read_user_input.h"
#include "headers/genRandMat.h"
#include "headers/eulerInt.h"
#include "headers/acorr_analysis.h"
#include "headers/other.h"

//std headers
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>
#include <random>

using namespace::std;
using namespace::Eigen;

int main(int argc, char* argv[])
{

	string name;	// output name

	int N;			// number of neurons
	double tf;		// end time of integration
	int tisave;		// number of timesteps to save
	double tinit;	// integrate until tinit befor "real" integration starts
	double dt;		// Stepsize for EulerIntergration

	int seedJ;		// seed for the random generator
	int seedI;
	double g;		// gain parameter

	// read valuese of the variables from input file
	read_input(N,g,tf,tisave,tinit,dt,seedJ,seedI,name, "input.txt");

	// defined in functions.h
	FI phi(1.,1);

	// output name
	if(name!="") name = "_"+name;

	// convert input to usefull variables
	double dtsave = tf/tisave;		//	time steps between saving
	double std = g/sqrt(N);

	// Random number generators
	// Ndist: standard normal distribution
	// Udist: Uniform distribution [0,1]
	default_random_engine generatorJ(seedJ);
	default_random_engine generatorI(seedI);
	normal_distribution<double> Ndist(0.0,1.0);
	uniform_real_distribution<double> Udist(0.0,1.0);


	// Initialize connectivity matrix J
	// and vector with E/I labels
	MatrixXd J(N,N);
	genRandMat(J,N,std,Ndist,generatorJ);

	// initialize x and r
	// x: synaptic potentials
	// r: firing rates
	VectorXd x(N);
	MatrixXd xt(N,tisave);
	VectorXd r(N);
//	MatrixXd rt(N,tisave);
	for(int i=0;i<N;++i) x(i) = (1./N)*Ndist(generatorI);
	phi(r,x);

	for(double t=0;t<tinit;t+=dt)
		eulerInt(x,r,J,phi,dtsave,dt);

	// start integrating
	for(int ti=0;ti<tisave;++ti) {
		eulerInt(x,r,J,phi,dtsave,dt);
		xt.col(ti) = x;
//		rt.col(ti) = r;
	}

/*
	End of the integration of the equations.
	From here only analysis.
*/

	// pow2 returns the smalles number that is a power 
	// of 2 and is larger than tisave. Muliplied by 2
	// for zero-padding for autocorrelations and PSD
	int t2 = 2*pow2(tisave);
	VecDoub delta(tisave,0.0);
	VecDoub psd_delta(t2/2,0.0);
//	VecDoub C(tisave,0.0);
//	VecDoub psd_C(t2/2,0.0);

//	acorr_analysis(xt,rt,delta,psd_delta,C,psd_C,N,tisave,t2);
	acorr_analysis(xt,delta,psd_delta,N,tisave,t2);

	// write results
	write_matrix(delta,tisave,"delta"+name+".csv");
//	write_matrix(C,tisave,"C"+name+".csv");
	write_matrix(psd_delta,t2/2,"psd_delta"+name+".csv");
//	write_matrix(psd_C,t2/2,"psd_C"+name+".csv");

	// transform delta to delta/delta0
	for(int ti=(tisave-1);ti>=0;--ti)
		delta[ti] /= delta[0];
	write_matrix(delta,tisave,"deltaN"+name+".csv");

	// claculate q=1-delta/delta0
	for(int ti=0;ti<tisave;++ti)
		delta[ti] = 1-delta[ti];
	write_matrix(delta,tisave,"q"+name+".csv");

	VecDoub average(tisave,0.0);
	for(int ti=0;ti<tisave;++ti) {
		for(int i=0;i<N;++i) {
			average[ti] += xt(i,ti)/double(N);
		}
	}
	write_matrix(average,tisave,"average"+name+".csv");


	VecDoub freq(t2/2,0.0);
	psd_frequencies(t2,dtsave,freq);

	vector<double> tval(tisave,0.0);
	for(int ti=0;ti<tisave;++ti)
		tval[ti] = ti*dtsave;

	// write  the rest of the results
	write_matrix(xt,N,tisave,"x" + name + ".csv");
//	write_matrix(rt,N,tisave,"phix"+name+".csv");
	write_matrix(tval,tisave,"t"+name+".csv");
	write_matrix(freq,t2/2,"f"+name+".csv");



	
	return 0;
}


