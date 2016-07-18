#ifndef GUARD_derivatives_h
#define GUARD_derivatives_h

#include <math.h>
#include "functions.h"

struct NW {
	std::vector<std::vector<double> >* w;
	int N;

	FI f;
	NW(std::vector<std::vector<double> >* ww, int NN, FI ff)
			: w(ww), N(NN), f(ff) {}

	void operator() (const Doub t, VecDoub_I& x, VecDoub_O& dxdt)
	{
		for(int i=0;i<N;++i) {
			dxdt[i] = -1.*x[i];
			for(int j=0;j<N;++j) {
				dxdt[i] += (*w)[i][j]*f(x[j]);
			}
		}
	}
};

#endif






