// nr3 headers
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
//#include "nr_headers/erf.h"
#include "nr_headers/fourier.h"
#include "nr_headers/correl.h"
#include "nr_headers/spectrum.h"
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
#include "headers/read_user_input.h"
#include "headers/functions.h"
#include "headers/other.h"
#include "headers/gen_rand_mat.h"

//std headers
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>

using namespace::std;


double nof(double a) { return 1.;}

int main(int argc, char* argv[])
{

	string name;	// output name

	int N;			// number of neurons
	int p;			// 1= one pupulation, 2= exc. and inh. pop.
	int function;	// which function to use
	int tf;			// end time of integration
	int tsave;		// number of timesteps to save
	int seed;		// seed for the random generator
	double g;		// gain parameter
	
	double meanE;	// mean of the exc. pop
	double meanI;	// mean of the inh. pop

	double a;		// if p=2 varE=g^2/(Na), varI=g^2/N
					// if p=1 a is not used

	double mean_noise;	// mean of the input noise
	double var_noise;	// var of the input noise

	int db;			// if db=1 detailed balence contition satisfied

	// read valuese of the variables from input file
	read_input(N,p,g,tf,tsave,function,mean_noise,var_noise,meanE,meanI,a,db,seed, "input.txt");

	// integration settings
	double atol;
	double rtol;
	double h1;
	double hmin;

	// read integration variables
	read_integration_vars(atol,rtol,h1,hmin,"integration_variables.txt");

	// over ride variable if user input is supplied
	if(argc >1){
		read_user_input(argc,argv,N,p,g,seed,tf,
			tsave,function,mean_noise,var_noise,
			meanE,meanI,a,db,name);
	}

	meanE /=sqrt(N);
	meanI /=sqrt(N);
	mean_noise /= sqrt(N);
	double std_noise = sqrt(var_noise);
	double dt = tf/(double)tsave;
	double sqrt_dt = sqrt(dt);
	double std  = g/sqrt((double)N);
	double stdE = g/sqrt(N*a);
	double stdI = g/sqrt(N);

	Ran r(seed);	// random number generator

	// f is a pointer to the transfer function
	double (*f)(double);
	if(function == 1) f = th_linear;
	if(function == 2) f = tanh;
	if(function == 3) f = tanh01;
	if(function == 4) f = tanh001;

	// create connectivity matrix w
//	vector<double> temp(N,.0);
	vector<vector<double> > w(N,vector<double> (N,0.0));

	// state vector of the neurons
	VecDoub x(N,0.0);
	vector<double> tval(tsave,0.0);
	vector<vector<double> > xt(N,tval);

	// initialize  x w rand uniform(-1,1)
	for(int i=0;i<N;++i) x[i] = 2*(0.5-r.doub());
	//initialize w
	if(p==1) gen_rand_mat(w,N,std,r);
	if(p==2) gen_rand_mat(w,N,meanE,meanI,stdE,stdI,db,r);


	// start integration
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


	int t2 = 2*pow2(tsave);
	// correlations
	VecDoub temp(t2,0.0);
	VecDoub acorr_temp(t2,0.0);
	VecDoub acorr(t2,0.0);
	VecDoub psd(t2,0.0);

	VecDoub mean(t2,0.0);
	VecDoub acorr_mean(t2,0.0);
	VecDoub psd_mean(t2,0.0);

	double mean_avg =0;
	for(int i=0;i<N;++i) {
		double temp_avg = 0;
		for(int ti=0;ti<tsave;++ti) {
			temp[ti] = xt[i][ti];
			temp_avg += temp[ti];

			mean[ti] += xt[i][ti]/(double)N;
		}

		for(int ti=0;ti<tsave;++ti)
			temp[ti] -= temp_avg/(double)tsave;

		correl(temp,temp,acorr_temp);
		
		for(int ti=0;ti<t2;++ti)
			acorr[ti] += acorr_temp[ti]/(N*acorr_temp[0]);

//		// psd 
//		realft(temp,1);
//		psd[0] += temp[0]/(double)N;
//		psd[tsave/2] += temp[1]/(double)N;
//		for(int ti=1;ti<(tsave-1);++ti)
//			psd[ti] += sqrt(temp[2*ti]*temp[2*ti]+temp[2*ti+1]*temp[2*ti+1])/N;
//		
	}

	for(int ti=0;ti<tsave;++ti) mean_avg += mean[ti];
	for(int ti=0;ti<tsave;++ti)
		mean[ti] -= mean_avg/(double)tsave;

	correl(mean,mean,acorr_mean);
	double acorr_mean0 = acorr_mean[0];
	for(int ti=0;ti<t2;++ti)
		acorr_mean[ti] /= acorr_mean0;

//	// psd_mean
//	realft(mean,1);
//	psd_mean[0] = mean[0];
//	psd_mean[tsave/2] = mean[1];
//	for(int ti=0;ti<(tsave-1);++ti)
//		psd_mean[ti] = sqrt(mean[2*ti]*mean[2*ti] + mean[2*ti+1]*mean[2*ti+1]);
//
//


	if(name!="") name = "_"+name;
//	write_matrix(psd,tsave/2,"psd"+name+".csv");
//	write_matrix(psd_mean,tsave/2,"psd_mean"+name+".csv");
	write_matrix(acorr,t2/2,"acorr" + name + ".csv");
	write_matrix(acorr_mean,t2/2,"acorr_mean"+name+".csv");
	write_matrix(tval,tsave,"t" + name + ".csv");
	write_matrix(xt,N,tsave,"x" + name + ".csv");
	write_matrix(mean,tsave,"mean"+name+".csv");	
	return 0;
}


