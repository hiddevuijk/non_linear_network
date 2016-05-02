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
	int p;			// 1= one pupulation, 2= exc. and inh. pop.
	double tf;		// end time of integration
	int tisave;		// number of timesteps to save
	double tinit;	// integrate until tinit befor "real" integration starts
	double dt;		// Stepsize for EulerIntergration

	int seed;		// seed for the random generator
	double g;		// gain parameter
	double sp; 		// sparseness (sp=0.1 -> 10% connections)
	double meanE;	// mean of the exc. pop
	double meanI;	// mean of the inh. pop
	double a;		// if p=2 varE=g^2/(Na), varI=g^2/N
					// if p=1 a is not used
	int db;			// if db=1 detailed balence contition satisfied
	double r0;		// variable of the transfer function; see functions.h.

	// read valuese of the variables from input file
	read_input(N,p,g,sp,tf,tisave,tinit,dt,r0,meanE,meanI,a,db,seed,name, "input.txt");


	// over ride variable if user input is supplied
	if(argc >1){
		read_user_input(argc,argv,N,p,g,sp,seed,tf,
			tisave,tinit,dt,r0,meanE,meanI,a,db,name);
	}

	// initialize transfer function functor
	// defined in functions.h
	FI phi(r0,p);

	// output name
	if(name!="") name = "_"+name;

	// convert input to usefull variables
	double dtsave = tf/tisave;		//	time steps between saving
	double f = meanI/(meanI-meanE);	// fraction of exc. pop.
	meanE /=sqrt(N);				// renormalized meanE
	meanI /=sqrt(N);				// renormalized meanI
	double stdE = g/sqrt(N*a);		// renormalized std of exc. pop.
	double stdI = g/sqrt(N);		// renormalized std of inh. pop.
	int Ne = f*N;					// number in exc. population
	int Ni =N-Ne; 					// number in inh. pop.


	// Random number generators
	// Ndist: standard normal distribution
	// Udist: Uniform distribution [0,1]
	default_random_engine generator(seed);
	normal_distribution<double> Ndist(0.0,1.0);
	uniform_real_distribution<double> Udist(0.0,1.0);


	// Initialize connectivity matrix J
	// and vector with E/I labels
	MatrixXd J(N,N);
	vector<char> EI(N,'E');
	if(p==2) {
		for(int i=Ne;i<N;++i)
			EI[i] = 'I';
		genRandMat(J,EI,N,meanE,meanI,stdE,stdI,db,Ndist,generator);
	} else {
		genRandMat(J,N,stdE,Ndist,generator);
	}

	// initialize x and r
	// x: synaptic potentials
	// r: firing rates
	VectorXd x(N);
	MatrixXd xt(N,tisave);
	VectorXd r(N);
	MatrixXd rt(N,tisave);
	for(int i=0;i<N;++i) x(i) = (1./N)*Ndist(generator);
	phi(r,x);

	// start integrating until tinit
	// dont save results in xt
	for(double t=0;t<tinit;t+=dt)
		eulerInt(x,r,J,phi,dtsave,dt);

	// start integrating
	for(int ti=0;ti<tisave;++ti) {
		eulerInt(x,r,J,phi,dtsave,dt);
		xt.col(ti) = x;
	}

/*
	End of the integration of the equations.
	From here only analysis.
*/

	// pow2 returns the smalles number that is a power 
	// of 2 and is larger than tisave. Muliplied by 2
	// for zero-padding for autocorrelations and PSD
	int t2 = 2*pow2(tisave);

	if(p==2) {
		VecDoub deltaE(tisave,0.0);
		VecDoub	deltaI(tisave,0.0);
		VecDoub psd_deltaE(t2/2,0.0);
		VecDoub psd_deltaI(t2/2,0.0);
		VecDoub CE(tisave,0.0);
		VecDoub CI(tisave,0.0);
		VecDoub psd_CE(t2/2,0.0);
		VecDoub psd_CI(t2/2,0.0);

		acorr_analysis(xt,rt,EI,deltaE,deltaI,psd_deltaE,psd_deltaI,CE,CI,psd_CE,psd_CI,N,Ne,Ni,tisave,t2);
		
		// write results
		write_matrix(deltaE,tisave,"deltaE"+name+".csv");
		write_matrix(deltaI,tisave,"deltaI"+name+".csv");
		write_matrix(CE,tisave,"CE"+name+".csv");
		write_matrix(CI,tisave,"CI"+name+".csv");
		write_matrix(psd_deltaE,t2/2,"psd_deltaE"+name+".csv");
		write_matrix(psd_deltaI,t2/2,"psd_deltaI"+name+".csv");
		write_matrix(psd_CE,t2/2,"psd_CE"+name+".csv");
		write_matrix(psd_CI,t2/2,"psd_CI"+name+".csv");

		// transform delta to delta/delta0
		for(int ti=(tisave-1);ti>=0;--ti){
			deltaE[ti] /= deltaE[0];
			deltaI[ti] /= deltaI[0];
		}
		write_matrix(deltaE,tisave,"deltaEN"+name+".csv");
		write_matrix(deltaI,tisave,"deltaIN"+name+".csv");
	
		// claculate q=1-delta/delta0
		for(int ti=0;ti<tisave;++ti){
			deltaE[ti] = 1-deltaE[ti];
			deltaI[ti] = 1-deltaI[ti];
		}
		write_matrix(deltaE,tisave,"qE"+name+".csv");
		write_matrix(deltaI,tisave,"qI"+name+".csv");

		VecDoub averageE(tisave,0.0);
		VecDoub averageI(tisave,0.0);
		for(int ti=0;ti<tisave;++ti) {
			for(int i=0;i<N;++i) {
				if(EI[i] == 'E')
					averageE[ti] += xt(i,ti)/double(Ne);
				if(EI[i] == 'I')
					averageI[ti] += xt(i,ti)/double(Ni);
			}
		}
		write_matrix(averageE,tisave,"averageE"+name+".csv");
		write_matrix(averageI,tisave,"averageI"+name+".csv");

	} else {
		VecDoub delta(tisave,0.0);
		VecDoub psd_delta(t2/2,0.0);
		VecDoub C(tisave,0.0);
		VecDoub psd_C(t2/2,0.0);

		acorr_analysis(xt,rt,delta,psd_delta,C,psd_C,N,tisave,t2);

		// write results
		write_matrix(delta,tisave,"delta"+name+".csv");
		write_matrix(C,tisave,"C"+name+".csv");
		write_matrix(psd_delta,t2/2,"psd_delta"+name+".csv");
		write_matrix(psd_C,t2/2,"psd_C"+name+".csv");

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
	}


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


