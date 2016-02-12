// nr3 headers
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
#include "nr_headers/fourier.h"
#include "nr_headers/ludcmp.h"
#include "nr_headers/correl.h"
#include "nr_headers/stepper.h"
#include "nr_headers/stepperbs.h"
#include "nr_headers/odeint.h"
#include "nr_headers/stepperdopr853.h"
#include "nr_headers/stepperdopr5.h"
#include "nr_headers/stepperstoerm.h"
#include "nr_headers/stepperross.h"
#include "nr_headers/steppersie.h"

// other headers
#include "headers/derivatives.h"
#include "headers/write_matrix.h"
#include "headers/box_muller.h"
#include "headers/read_input.h"
#include "headers/read_user_input.h"
#include "headers/functions.h"
#include "headers/other.h"
#include "headers/gen_rand_mat.h"
#include "headers/acorr.h"
#include "headers/gc.h"

//std headers
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>

using namespace::std;

int main(int argc, char* argv[])
{

	string name;	// output name

	int N;			// number of neurons
	int p;			// 1= one pupulation, 2= exc. and inh. pop.
	int tf;			// end time of integration
	int tsave;		// number of timesteps to save
	int tinit;		// integrate until tinit befor "real" integration starts

	int seed;		// seed for the random generator
	double g;		// gain parameter
	
	double meanE;	// mean of the exc. pop
	double meanI;	// mean of the inh. pop

	double a;		// if p=2 varE=g^2/(Na), varI=g^2/N
					// if p=1 a is not used


	double mean_noise;	// mean of the input noise
	double var_noise;	// var of the input noise

	int db;			// if db=1 detailed balence contition satisfied

	double r0;		// r0 value in F-I functor ( def in function.h)

	// read valuese of the variables from input file
	read_input(N,p,g,tf,tsave,tinit,r0,mean_noise,var_noise,meanE,meanI,a,db,seed,name, "input.txt");

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
			tsave,tinit,r0,mean_noise,var_noise,
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

	Tanhr0 f(r0,p);

	// create connectivity matrix w
	vector<vector<double> > w(N,vector<double> (N,0.0));
	vector<vector<double> >* wptr = &w;
	vector<char> EI(N,'E');

	// state vector of the neurons
	VecDoub x(N,0.0);
	vector<double> tval(tsave,0.0);
	vector<vector<double> > xt(N,tval);
	vector<vector<double> > xtphi(N,tval);
	vector<vector<double> > dh(N,tval);

	// initialize  x w rand uniform(-1,1)
	for(int i=0;i<N;++i) x[i] = 2*(0.5-r.doub());
	//initialize w
	if(p==2) gen_rand_mat(w,N,meanE,meanI,stdE,stdI,db,EI,r);
	else gen_rand_mat(w,N,std,r);

	// start integration
	NW nw(wptr,N,f);
//	NW_ross nwr(w,N,f);
	Output out;

	Output out_init;
	Odeint<StepperDopr853<NW> > ode_start(x,0,tinit,atol,rtol,h1,hmin,out_init,nw);

	for(int ti=0;(ti+1)<tsave;++ti){
		Odeint<StepperDopr853<NW> > ode(x,ti*dt,(ti+1)*dt,atol,rtol,h1,hmin,out,nw);
//		Odeint<StepperRoss<NW_ross> > ode(x,ti*dt,(ti+1)*dt,atol,rtol,h1,hmin,out,nwr);
		ode.integrate();
		for(int i=0;i<N;++i) {
			xt[i][ti] = x[i];
			xtphi[i][ti] = f(xt[i][ti]);

			x[i] += sqrt_dt*std_noise*bm_transform(r);
		}
		tval[ti]= ti*dt;
	}
	tval[tsave-1] = tf;


/*
	End of the integration of the equations.
	From here only analysis.
*/

	// pow2 returns the smalles number that is a power 
	// of 2 and is larger than tsave. Muliplied by 2
	// for zero-padding for autocorrelations and PSD
	int t2 = 2*pow2(tsave);

	// correlations
	VecDoub acorrn(tsave,0.0);
	VecDoub	acorr(t2,0.0);
	VecDoub psd(t2/2,0.0);
	VecDoub acorrn_phi(tsave,0.0);
	VecDoub acorr_phi(tsave,0.0);
	VecDoub psd_phi(tsave,0.0);
	VecDoub freq(t2/2,0.0);

	VecDoub dh_temp(t2,0.0);
	VecDoub acorr_dh_temp(t2,0.0);
	VecDoub acorr_dh(tsave,0.0);
	VecDoub psd_dh(tsave,0.0);
	VecDoub psd_dh_temp(t2,0.0);

	// first mean then ...
	VecDoub average_x(t2,0.0);
	VecDoub acorrn_mean(t2,0.0);
	VecDoub acorr_mean(t2,0.0);
	VecDoub psd_mean(t2/2,0.0);
	VecDoub average_phi(t2,0.0);
	VecDoub acorrn_mean_phi(t2,0.0);
	VecDoub acorr_mean_phi(t2,0.0);
	VecDoub psd_mean_phi(t2/2,0.0);

	// temporary vec
	VecDoub temp(t2,0.0);
	VecDoub temp_phi(t2,0.0);

	VecDoub acorrn_temp(t2,0.0);
	VecDoub acorr_temp(t2,0.0);
	VecDoub psd_temp(t2/2,0.0);
	VecDoub acorrn_temp_phi(t2,0.0);
	VecDoub acorr_temp_phi(t2,0.0);
	VecDoub psd_temp_phi(t2,0.0);

	//population means
	double pmE = 0;
	double pmI = 0;	
	for(int i=0;i<N;++i){
		for(int ti=0;ti<tsave;++ti) {
			if(EI[i]=='E') pmE += xt[i][ti]/(N*tsave);
			else pmI += xt[i][ti]/(N*tsave);
		}
	}

	// pmi is either pmE or pmI, depending on i
	double pmi = 0;

	for(int i=0;i<N;++i) {
		if(EI[i] == 'E') pmi = pmE;
		else pmi = pmI;

		for(int ti=0;ti<tsave;++ti){
			temp[ti] = xt[i][ti];
			temp_phi[ti] = xtphi[i][ti];

			average_x[ti] += temp[ti]/(double)N;		
			average_phi[ti] += temp_phi[ti]/(double)N;

			dh_temp[ti] = temp[ti] - pmi;
			dh[i][ti] = dh_temp[ti];
		}

		autocorrel_norm(temp,tsave,acorrn_temp);
		autocorrel_psd(temp,tsave,acorr_temp,psd_temp);

		autocorrel_norm(temp_phi,tsave,acorrn_temp_phi);
		autocorrel_psd(temp_phi,tsave,acorr_temp_phi,psd_temp_phi);

		autocorrel_psd(dh_temp,tsave,acorr_dh_temp,psd_dh_temp);
		double dh_temp0 = dh_temp[0];		

		for(int ti=0;ti<tsave;++ti){
			acorrn[ti] += acorrn_temp[ti]/(double)N;
			acorr[ti] += acorr_temp[ti]/(double)N;

			acorrn_phi[ti] += acorrn_temp_phi[ti]/(double)N;
			acorr_phi[ti] += acorr_temp_phi[ti]/(double)N;

			acorr_dh[ti] += acorr_dh_temp[ti]/(dh_temp0*N);

			if(ti<t2/2)	{
				psd[ti] += psd_temp[ti]/(double)N;
				psd_phi[ti] += psd_temp_phi[ti]/(double)N;
				psd_dh[ti] += psd_dh_temp[ti]/(double)N;
			}
		}

	}	

	autocorrel_norm(average_x,tsave,acorrn_mean);
	autocorrel_psd(average_x,tsave,acorr_mean,psd_mean);

	autocorrel_norm(average_phi,tsave,acorrn_mean_phi);
	autocorrel_psd(average_phi,tsave,acorr_mean_phi,psd_mean_phi);

	psd_frequencies(t2,dt,freq);

	// write results
	if(name!="") name = "_"+name;
	write_matrix(xt,N,tsave,"x" + name + ".csv");
	write_matrix(xtphi,N,tsave,"phix"+name+".csv");
	write_matrix(tval,tsave,"t"+name+".csv");
	write_matrix(average_x,tsave,"mean"+name+".csv");	
	write_matrix(average_phi,tsave,"mean_phi"+name+".csv");
	write_matrix(freq,t2/2,"freq"+name+".csv");
	write_matrix(psd,t2/2,"psd"+name+".csv");
	write_matrix(psd_phi,t2/2,"psd_phi"+name+".csv");
	write_matrix(psd_mean,t2/2,"psd_mean"+name+".csv");
	write_matrix(psd_mean_phi,t2/2,"psd_mean_phi"+name+".csv");
	write_matrix(acorrn,tsave,"acorrn"+name+".csv");
	write_matrix(acorrn_phi,tsave,"acorrn_phi"+name+".csv");
	write_matrix(acorrn_mean,tsave,"acorrn_mean"+name+".csv");
	write_matrix(acorrn_mean_phi,tsave,"acorrn_mean_phi"+name+".csv");
	write_matrix(acorr,tsave,"acorr"+name+".csv");
	write_matrix(acorr_phi,tsave,"acorr_phi"+name+".csv");
	write_matrix(acorr_mean,tsave,"acorr_mean"+name+".csv");
	write_matrix(acorr_mean_phi,tsave,"acorr_mean_phi"+name+".csv");
	write_matrix(acorr_dh,tsave,"acorr_dh"+name+".csv");
	write_matrix(psd_dh,tsave,"psd_dh"+name+".csv");
	write_matrix(dh,N,tsave,"dh"+name+".csv");	
	return 0;
}
