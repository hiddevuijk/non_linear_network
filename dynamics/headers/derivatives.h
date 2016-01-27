#ifndef GUARD_derivatives_h
#define GUARD_derivatives_h

#include <math.h>


struct NW {
	vector<vector<double> > w;
	int N;
	double (*f)(double);
	NW(vector<vector<double> > ww,int  NN, double (*ff)(double)) :
		w(ww) , N(NN),f(ff) {}

	void operator() (const Doub t, VecDoub_I& x, VecDoub_O& dxdt)
	{
		for(int i=0;i<N;++i) {
			dxdt[i] = -1.*x[i];
			for(int j=0;j<N;++j) {
				dxdt[i] += w[i][j]*f(x[j]);
			}
		}
	}
};



#endif
