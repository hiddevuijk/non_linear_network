#ifndef GUARD_functions_h
#define GUARD_functions_h

#include <string>
#include <math.h>

// transfer functions

struct FI {
private:
	double v;
public:
	FI(double vv) : v(vv) {}

	double operator()(double x){
		if(x<=0) return 0.;
		 else return pow(x,v);
	}
};





#endif
