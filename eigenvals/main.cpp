#include </usr/local/Cellar/eigen/3.2.4/include/eigen3/Eigen/Dense>

#include "nr_headers/nr3.h"
#include "nr_headers/ran.h"
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
	double var; 	// variance units: 1/N

	int seed;		//seed for the random number generator
	
	read_input(N,meanE,meanI,var,seed, "input.txt");

	Ran r(seed);

	// prob. of exc. col. in w
	// s.t. f*meanE + (1-f)*meanI =0
	double f;
	if( meanI-meanE == 0) f = 0.5;
	else f = meanI/(meanI-meanE);

	double std = sqrt(var/(double)N);
	meanE = meanE/sqrt((double)N);
	meanI = meanI/sqrt((double)N);

	vector<double> m(N,meanI);
	for(int i =0;i<f*N;++i) m[i] = meanE;
	shuffle(m,N,r);
	double sum =0;
	for(int i=0;i<N;++i) sum+=m[i];
	cout << sum << endl;


	// create w, the matrix of connection strengths
	MatrixXcf w(N,N);
	vector<double> row_sum(N,0.0);
	for(int i=0;i<N;++i){
		for( int j=0;j<N;++j) {
			w(i,j) = std*bm_transform(r);
			row_sum[i] += w(i,j).real();
		}
	}
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			w(i,j) += m[j]-1.*row_sum[i]/N;
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


