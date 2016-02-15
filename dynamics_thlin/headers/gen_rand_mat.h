#ifndef GUARD_gen_rand_mat_h
#define GUARD_gen_rand_mat_h

#include <vector>
#include "box_muller.h"
#include "shuffle.h"

void gen_rand_mat(std::vector<std::vector<double> >& w,
		int N, double std,double gm, Ran& r)
{
	for(int i=0;i<N;++i) {
		for(int j=0;j<N;++j) {
			w[i][j] = gm+std*bm_transform(r);
		}
	}
}



#endif
