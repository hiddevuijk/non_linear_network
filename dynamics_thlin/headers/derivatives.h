#ifndef GUARD_derivatives_h
#define GUARD_derivatives_h

#include <math.h>
#include "functions.h"

struct NW {
	std::vector<std::vector<double> >* w;
	int N;
	double h0;
	FI f;
	NW(std::vector<std::vector<double> >* ww, int NN, FI ff, double h00)
			: w(ww), N(NN), f(ff), h0(h00) {}

	void operator() (const Doub t, VecDoub_I& x, VecDoub_O& dxdt)
	{
		for(int i=0;i<N;++i) {
			dxdt[i] = -1.*x[i]+h0;
			for(int j=0;j<N;++j) {
				dxdt[i] += (*w)[i][j]*f(x[j]);
			}
		}
	}
};

#endif






