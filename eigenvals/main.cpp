#include </usr/local/Cellar/eigen/3.2.4/include/eigen3/Eigen/Dense>

#include <random>
#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
#include "nr_headers/erf.h"
#include "nr_headers/gamma.h"
#include "nr_headers/deviates.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>

using namespace::std;
using namespace::Eigen;

#include "headers/write_matrix.h"
#include "headers/box_muller.h"
#include "headers/read_input.h"
#include "headers/shuffle.h"

int main()
{

	int N;			// number of neurons
	double meanE;	// mean strength of exc. neurons in units 1/sqrt(N)
	double meanI;	// mean strength of inh. neurons in units 1/sqrt(N)
	double g; 		// gain parameter
	double a;		// varE= g^2/(Na), varI=g^2/N
	int db;
	int seed;		//seed for the random number generator

	string dist;	// normal or  lognormal
		
	read_input(N,meanE,meanI,g,a,db,dist,seed, "input.txt");

	Ran r(seed);
	default_random_engine generator;
	

	// prob. of exc. col. in w
	// s.t. f*meanE + (1-f)*meanI =0
	double f;
	if( meanI-meanE == 0) f = 0.5;
	else f = meanI/(meanI-meanE);


	double varE = g*g/(a*N);
	double varI = g*g/N;
	double stdE = g/sqrt(N*a);
	double stdI = g/sqrt((double)N);
	meanE = meanE/sqrt((double)N);
	meanI = meanI/sqrt((double)N);

	normal_distribution<double> normdist01(0.0,1.0);
	normal_distribution<double> normdistE(meanE,stdE);
	normal_distribution<double> normdistI(meanI,stdI);

	double muE = log(meanE/sqrt(1+varE/(meanE*meanE)));
	double muI = log(-1.*meanI/sqrt(1+varI/(meanI*meanI)));
	double sigE = sqrt(log(1+varE/(meanE*meanE)));
	double sigI = sqrt(log(1+varI/(meanI*meanI)));


	lognormal_distribution<double> lognormdistE(muE,sigE);
	lognormal_distribution<double> lognormdistI(muI,sigI);



	vector<char> EI(N,'I');
	for(int i=0;i<f*N;++i) EI[i] = 'E';
	shuffle(EI,N,r);


	// create w, the matrix of connection strengths
	MatrixXcf w(N,N);
	vector<double> row_sum(N,0.0);


	if(dist=="normal") {	
		// row=i, col=j
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				if(EI[j]=='E') {
					w(i,j) = stdE*normdist01(generator);
					row_sum[i] += w(i,j).real();
				}
				if(EI[j]=='I') {
					w(i,j) = stdI*normdist01(generator);
					row_sum[i] += w(i,j).real();
				}
			}
		}
		// if db = 1 impose detailed balance condition
		if(db==1) {
			for(int i=0;i<N;++i) {
				for(int j=0;j<N;++j) {
					w(i,j) -= row_sum[i]/N;
				}
			}
		}

		// add M to W
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				if(EI[j]=='E') w(i,j) +=  meanE;
				if(EI[j]=='I') w(i,j) +=  meanI;
			}
		}
	}

	if(dist=="lognormal") {	
		// row=i, col=j
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				if(EI[j]=='E') {
					w(i,j) = lognormdistE(generator);
				}
				if(EI[j]=='I') {
					w(i,j) = -1.*lognormdistI(generator);
				}
			}
		}
	}

	// Eigen object for computing eigensystem
	ComplexEigenSolver<MatrixXcf> eigenw;
	eigenw.compute(w);

	vector<double> eval_re(N);
	vector<double> eval_im(N);

	for(int i=0;i<N;i++) {
		eval_re[i] = eigenw.eigenvalues()[i].real();
		eval_im[i] = eigenw.eigenvalues()[i].imag();
	}

	write_matrix(eval_re,N,"eval_re.csv");	
	write_matrix(eval_im,N,"eval_im.csv");

	return 0;
}


