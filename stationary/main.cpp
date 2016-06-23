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
#include "headers/getStat.h"
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

	double tf;		// end time of integration
	int tisave;		// number of timesteps to save
	int Nstat;		// number of time to in tegrate tf and calculate stats
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
	read_input(N,p,g,tf,tisave,Nstat,dt,r0,meanE,meanI,a,db,seed,name, "input.txt");

	if(name!="") name = "_"+name;

	double dtsave = tf/(double)tisave;
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
		write_matrix(EI,EI.size(),"EI"+name+".csv");
	}
	//initialize x and r
	VectorXd x(N);
	VectorXd r(N);
	for(int i=0;i<N;++i) x(i) = (1./N)*Ndist(generator);
	phi(r,x);


	MatrixXd m0(N,Nstat);
	MatrixXd m1(N,Nstat);
	MatrixXd m2(N,Nstat);
	
	if(p==2) {
		MatrixXd xEt(Ne,tisave);
		MatrixXd xIt(Ni,tisave);
		for(int Nstati=0;Nstati<Nstat;++Nstati) {
			for(int ti=0;ti<tisave;++ti) {
				eulerInt(x,r,J,phi,dtsave,dt);
				xEt.col(ti) = x.head(Ne);
				xIt.col(ti) = x.tail(Ni);
			}
			getStat(xEt,xIt,m0,m1,m2,Nstati);	
		}	
		//write m0,m1,m2
		write_matrix(m0,m0.rows(),m0.cols(),"m0.csv");
		write_matrix(m1,m1.rows(),m1.cols(),"m1.csv");
		write_matrix(m2,m2.rows(),m2.cols(),"m2.csv");
	} else {
		MatrixXd xt(N,tisave);
		for(int Nstati=0;Nstati<Nstat;++Nstati) {
			for(int ti=0;ti<tisave;++ti) {
				eulerInt(x,r,J,phi,dtsave,dt);
				xt.col(ti) = x;
			}
			getStat(xt,m0,m1,m2,Nstati);
		}
		//write m0,m1,m2
		write_matrix(m0,m0.rows(),m0.cols(),"m0.csv");	
		write_matrix(m1,m1.rows(),m1.cols(),"m1.csv");
		write_matrix(m2,m2.rows(),m2.cols(),"m2.csv");
	}
	return 0;
}


