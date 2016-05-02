#include "Eigen/Core"

// other headers
#include "headers/write_matrix.h"
#include "headers/read_input.h"
#include "headers/functions.h"
#include "headers/other.h"
#include "headers/genRandMat.h"
#include "headers/eulerInt.h"

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
	int Ne = N; // number in exc. population
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
	write_matrix(EI,EI.size(),"EI_"+name+".csv");
	//initialize x and r
	VectorXd x(N);
	VectorXd r(N);
	for(int i=0;i<N;++i) x(i) = (1./N)*Ndist(generator);
	phi(r,x);

	// block for initial integration
	// save results
	{
		MatrixXd xtInit(N,tiInit);
		for(double ti=0;ti<tiInit;++ti) {
			eulerInt(x,r,J,phi,dtsave,dt);
			xtInit.col(ti) = x;
		}
		write_matrix(xtInit,N,tiInit,"xinit_"+name+".csv");
	}
	{
		MatrixXd xt(N,tisave);
		MatrixXd rt(N,tisave);

		for(int ti=0;ti<tisave;++ti) {
			eulerInt(x,r,J,phi,dtsave,dt);
			xt.col(ti) = x;
	//		rt.col(ti) = r;
		}

		// write  the rest of the results
		write_matrix(xt,N,tisave,"x_" + name + ".csv");
	}
	return 0;
}


