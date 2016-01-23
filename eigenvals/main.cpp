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

int main()
{

	int N;
	double meanS,meanA;
	double varS,varA;
	int seed;
	
	read_input(N,meanS,varS,meanA,varA,seed, "input.txt");
	Ran r(seed);

	MatrixXcf w(N,N);
	double mean =0;
	for(int i=0;i<N;++i){
		for( int j=0;j<i;++j) {
			double a = meanA+varA*bm_transform(r);
			double s = meanS+varS*bm_transform(r);
			w(i,j) = a+s;
			w(j,i) = -1.*a+s;
			mean += 2*s;
		}
		double s= meanS+varS*bm_transform(r);
		w(i,i) = s;
		mean += s;
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


