#ifndef GUARD_derivatives_h
#define GUARD_derivatives_h

#include <math.h>
#include "functions.h"

struct NW {
	std::vector<std::vector<double> >* w;
	int N;

	Tanhr0 f;
	NW(std::vector<std::vector<double> >* ww, int NN, Tanhr0 ff)
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

/// checl before use!!!!!
struct NW_ross {
	std::vector<std::vector<double> > w;
	int N;

	Tanhr0 f;
	NW_ross(std::vector<std::vector<double> > ww, int NN, Tanhr0 ff)
			: w(ww), N(NN), f(ff) {}

	void operator() (const Doub t, VecDoub_I& x, VecDoub_O& dxdt)
	{
		for(int i=0;i<N;++i) {
			dxdt[i] = -1.*x[i];
			for(int j=0;j<N;++j) {
				dxdt[i] += w[i][j]*f(x[j]);
			}
		}
	}

	void jacobian(const Doub t, VecDoub_I& x, VecDoub_O& dfdt,
			MatDoub_O& dfdx) {
		Int n=x.size();
		for(int i=0;i<n;++i){
			dfdt[i] =0.0;
			for(int j=0;j<n;++j)
				dfdx[i][j] = w[i][j] * f.derivative(x[j]);
			dfdx[i][i] -= 1.;
		}
	}

};

#endif






