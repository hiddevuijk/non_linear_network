#ifndef GUARD_gen_rand_mat_h
#define GUARD_gen_rand_mat_h

#include <vector>
#include "box_muller.h"

void gen_rand_mat(std::vector<std::vector<double> >& w,
		int N, double std, Ran r)
{
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			w[i][j] = std*bm_transform(r);
		}
	}
}

//void gen_rand_mat(std::vector<std::vector<double> >& w,
//		int N, double mean,double meanE, double meanI,
//		double stdE,double stdI,r)
//
//}
//

#endif
