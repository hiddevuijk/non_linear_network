#ifndef GUARD_gen_rand_mat_h
#define GUARD_gen_rand_mat_h

#include <vector>
#include "box_muller.h"
#include "shuffle.h"

void gen_rand_mat(std::vector<std::vector<double> >& w,
		int N, double std, Ran& r)
{
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			w[i][j] = std*bm_transform(r);
		}
	}
}

void gen_rand_mat(std::vector<std::vector<double> >& w,
				int N, double meanE, double meanI, 
				double stdE, double stdI, int db,Ran& r)
{
	double f;
	if( meanI-meanE == 0) f = 0.5;
	else f = meanI/(meanI-meanE);

	vector<char> EI(N,'I');
	for(int i=0;i<f*N;++i) EI[i] = 'E';
	shuffle(EI,N,r);


	// create w, the matrix of connection strengths
	std::vector<double> row_sum(N,0.0);
	
	// row=i, col=j
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			if(EI[j]=='E') {
				w[i][j] = stdE*bm_transform(r);
				row_sum[i] += w[i][j];
			}
			if(EI[j]=='I') {
				w[i][j] = stdI*bm_transform(r);
				row_sum[i] += w[i][j];
			}
		}
	}
	// if db = 1 impose detailed balance condition
	if(db==1) {
		for(int i=0;i<N;++i) {
			for(int j=0;j<N;++j) {
				w[i][j] -= row_sum[i]/N;
			}
		}
	}

	// add M to W
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			if(EI[j]=='E') w[i][j] += meanE;
			if(EI[j]=='I') w[i][j] += meanI;
		}
	}

}


#endif
