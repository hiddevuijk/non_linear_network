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

	int N=1000;			// number of neurons
	double g=1; 		// gain parameter
	int seed;		//seed for the random number generator

		

	Ran r(seed);
	default_random_engine generator;
	


	double std = g/sqrt(N);
	
	normal_distribution<double> normDist(0.0,1.0);




	// create w, the matrix of connection strengths
	MatrixXcf w(N,N);

	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			w(i,j) = std*normDist(generator);
		}
		w(i,i) -= 1;
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


