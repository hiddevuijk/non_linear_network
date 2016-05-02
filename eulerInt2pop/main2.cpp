#include "Eigen/Core"

// nr3 headers
#include "nr_headers/nr3.h"
#include "nr_headers/fourier.h"

// other headers
#include "headers/write_matrix.h"
#include "headers/read_input.h"
#include "headers/functions.h"
#include "headers/other.h"
#include "headers/genRandMat.h"
#include "headers/eulerInt.h"
#include "headers/acorr.h"
#include "headers/acorr_analysis.h"
#include "headers/mean.h"
#include "headers/other.h"
#include "headers/normalize.h"
//std headers
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>
#include <random>

// nr3 headers
#include "nr_headers/nr3.h"


using namespace::std;
using namespace::Eigen;

int main(int argc, char* argv[])
{

	string name;	// output name

	int N;			// number of neurons
	int p;			// 1= one pupulation, 2= exc. and inh. pop.
	double tf;			// end time of integration
	int tisave;		// number of timesteps to save
	double tinit;	// integrate until tinit befor "real" integration starts
	double dt;		// integration stepsize
	int seed;		// seed for the random generator
	double g;		// gain parameter
	double meanE;	// mean of the exc. pop
	double meanI;	// mean of the inh. pop
	double a;		// if p=2 varE=g^2/(Na), varI=g^2/N
					// if p=1 a is not used
	int db;			// if db=1 detailed balence contition satisfied
	double r0;		// parameter in the firing rate function

	// read valuese of the variables from input file
	read_input(N,p,g,tf,tisave,tinit,dt,r0,meanE,meanI,a,db,seed,name, "input.txt");

	if(name!="") name = "_"+name;

	double dtsave = tf/(double)tisave;
	int tiInit = tinit/dtsave;
	meanE /=sqrt(N);
	meanI /=sqrt(N);
	double std  = g/sqrt((double)N);
	double stdE = g/sqrt(N*a);
	double stdI = g/sqrt(N);
	int Ne = 0; // number in exc. population
	int Ni =0; // number in inh. pop.
	
	// random generator
	default_random_engine generator(seed);
	normal_distribution<double> Ndist(0.0,1.0);
	uniform_real_distribution<double> Udist(0.0,1.0);

	// firing rate function
	FI phi(r0,p);

	// exc. inh. neurons
	vector<char> EI(N,'E');

	// initialize J
	MatrixXd J(N,N);
	if(p==1) {
		genRandMat(J,N,std,Ndist,generator);
	} else {
		genRandMat(J,EI,N,meanE,meanI,stdE,stdI,db,Ndist,generator);
	}
	if(p==2) {	
		for(int i=0;i<N;++i) {
			if(EI[i]=='E')
				++Ne;
			else
				break;
		}
		Ni = N-Ne;
	}	

	write_matrix(EI,EI.size(),"EI"+name+".csv");
	//initialize x and r
	VectorXd x(N);
	VectorXd r(N);
	for(int i=0;i<N;++i) x(i) = (1./N)*Ndist(generator);
	phi(r,x);

	// block for initial integration
	// save results
	{
		MatrixXd xEtInit(Ne,tiInit);
		MatrixXd xItInit(Ni,tiInit);
		vector<double> tin(tiInit);
		for(int ti=0;ti<tiInit;++ti) {
			eulerInt(x,r,J,phi,dtsave,dt);
			xEtInit.col(ti) = x.head(Ne);
			xItInit.col(ti) = x.tail(Ni);
			tin[ti] = ti*dtsave;
		}
		write_matrix(xEtInit,Ne,tiInit,"xEinit"+name+".csv");
		write_matrix(xItInit,Ni,tiInit,"xIinit"+name+".csv");
		write_matrix(tin,tiInit,"t_init"+name+".csv");

		// analysis
		int t2 = 2*pow2(tiInit);
		
		VecDoub deltaE(tiInit,0.0);
		VecDoub deltaE_psd(t2/2,0.0);
		acorr_analysis(xEtInit,deltaE,deltaE_psd,Ne,tiInit,t2);
		write_matrix(deltaE,tiInit,"deltaE_init"+name+".csv");
		write_matrix(deltaE_psd,tiInit,"deltaE_psd_init"+name+".csv");
		normalize(deltaE,tiInit);
		write_matrix(deltaE,tiInit,"deltaEN_init"+name+".csv");


		VecDoub deltaI(tiInit,0.0);
		VecDoub deltaI_psd(t2/2,0.0);
		acorr_analysis(xItInit,deltaI,deltaI_psd,Ni,tiInit,t2);
		write_matrix(deltaI,tiInit,"deltaI_init"+name+".csv");
		write_matrix(deltaI_psd,tiInit,"deltaI_psd_init"+name+".csv");
		normalize(deltaI,tiInit);
		write_matrix(deltaI,tiInit,"deltaIN_init"+name+".csv");
		//averages
		VecDoub avgE(t2,0.0);
		mean(xEtInit,avgE);
		VecDoub avgI(t2,0.0);
		mean(xItInit,avgI);
		write_matrix(avgE,tiInit,"avgE_init"+name+".csv");
		write_matrix(avgI,tiInit,"avgI_init"+name+".csv");

		VecDoub acorrAvgE(t2,0.0);
		autocorrel_norm(avgE,Ne,acorrAvgE);
		write_matrix(acorrAvgE,tiInit,"acorrAvgE_init"+name+".csv");
		normalize(acorrAvgE,tiInit);
		write_matrix(acorrAvgE,tiInit,"acorrAvgEN_init"+name+".csv");

		VecDoub acorrAvgI(t2,0.0);
		autocorrel_norm(avgI,Ni,acorrAvgI);
		write_matrix(acorrAvgI,tiInit,"acorrAvgI_init"+name+".csv");
		normalize(acorrAvgI,tiInit);
		write_matrix(acorrAvgI,tiInit,"acorrAvgIN_init"+name+".csv");


	}
	{
		
		MatrixXd xEt(Ne,tisave);
		MatrixXd xIt(Ni,tisave);
		vector<double> t(tisave,0.0);
		for(int ti=0;ti<tisave;++ti) {
			eulerInt(x,r,J,phi,dtsave,dt);
			xEt.col(ti) = x.head(Ne);
			xIt.col(ti) = x.tail(Ni);
			t[ti] = ti*dtsave;
		}

		// write  the rest of the results
		write_matrix(xEt,Ne,tisave,"xE" + name + ".csv");
		write_matrix(xIt,Ni,tisave,"xI" + name + ".csv");
		write_matrix(t,tisave,'t'+name+".csv");

		// analysis
		int t2 = 2*pow2(tisave);
		
		VecDoub deltaE(tisave,0.0);
		VecDoub deltaE_psd(t2/2,0.0);

		acorr_analysis(xEt,deltaE,deltaE_psd,Ne,tisave,t2);
		write_matrix(deltaE,tisave,"deltaE"+name+".csv");
		write_matrix(deltaE_psd,tisave,"deltaE_psd"+name+".csv");
		normalize(deltaE,tisave);
		write_matrix(deltaE,tiInit,"deltaEN"+name+".csv");


		VecDoub deltaI(tisave,0.0);
		VecDoub deltaI_psd(t2/2,0.0);
		acorr_analysis(xIt,deltaI,deltaI_psd,Ni,tisave,t2);
		write_matrix(deltaI,tisave,"deltaI"+name+".csv");
		write_matrix(deltaI_psd,tisave,"deltaI_psd"+name+".csv");
		normalize(deltaI,tiInit);
		write_matrix(deltaI,tiInit,"deltaIN"+name+".csv");


		//averages
		VecDoub avgE(t2,0.0);
		mean(xEt,avgE);
		VecDoub avgI(t2,0.0);
		mean(xIt,avgI);
		write_matrix(avgE,tisave,"avgE"+name+".csv");
		write_matrix(avgI,tisave,"avgI"+name+".csv");
		VecDoub acorrAvgE(t2,0.0);
		autocorrel_norm(avgE,Ne,acorrAvgE);
		write_matrix(acorrAvgE,tisave,"acorrAvgE"+name+".csv");
		normalize(acorrAvgE,tiInit);
		write_matrix(acorrAvgE,tiInit,"acorrAvgEN"+name+".csv");

		VecDoub acorrAvgI(t2,0.0);
		autocorrel_norm(avgI,Ni,acorrAvgI);
		write_matrix(acorrAvgI,tisave,"acorrAvgI"+name+".csv");
		normalize(acorrAvgI,tiInit);
		write_matrix(acorrAvgI,tiInit,"acorrAvgIN"+name+".csv");

	}
	return 0;
}


