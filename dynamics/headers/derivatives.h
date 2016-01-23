#ifndef GUARD_derivatives_h
#define GUARD_derivatives_h

#include <math.h>


struct NW {
	vector<vector<double> > w;
	int N;
	NW(vector<vector<double> > ww,int  NN) :
		w(ww) , N(NN) {}

	void operator() (const Doub t, VecDoub_I& x, VecDoub_O& dxdt)
	{
		for(int i=0;i<N;++i) {
			dxdt[i] = -1.*x[i];
			for(int j=0;j<N;++j) {
				dxdt[i] += w[i][j]*tanh(x[j]);
			}
		}
	}
};



#endif
