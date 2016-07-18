#ifndef GUARD_genRandMat_h
#define GUARD_genRandMat_h

#include "../Eigen/Core"
#include <random>

void genRandMat(Eigen::MatrixXd& J, int N, double std,
		std::normal_distribution<double>& Ndist,
		std::default_random_engine& generator)
{
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			J(i,j) = std*Ndist(generator);
		}
	}

}

void genRandMat(Eigen::MatrixXd& J,std::vector<char>& EI,
				int N, double meanE, double meanI, double stdE,
				double stdI, int db,
				std::normal_distribution<double>& Ndist,
				std::default_random_engine& generator)
{

	double f;
	if(meanE-meanI == 0) f = 0.5;
	else f = meanI/(meanI-meanE);
	for(int i=0;i<f*N;++i) EI[i] = 'E';
	for(int i=f*N;i<N;++i) EI[i] = 'I';

	// create J, the matrix of connection strengths
	std::vector<double> row_sum(N,0.0);
	
	// row=i, col=j
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			if(EI[j]=='E') {
				J(i,j) = stdE*Ndist(generator);
				row_sum[i] += J(i,j);
			}
			if(EI[j]=='I') {
				J(i,j) = stdI*Ndist(generator);
				row_sum[i] += J(i,j);
			}
		}
	}
	// if db = 1 impose detailed balance condition
	if(db==1) {
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				J(i,j) -= row_sum[i]/N;
			}
		}
	}

	// add M to J
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			if(EI[j]=='E') J(i,j) += meanE;
			if(EI[j]=='I') J(i,j) += meanI;
		}
	}

}


#endif
