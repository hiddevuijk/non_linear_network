#ifndef GUARD_functions_h
#define GUARD_functions_h

#include <string>
#include <math.h>

// transfer functions

struct FI {
private:
	double v;
	std::string f;
public:
	FI(double vv, std::string ff) : v(vv), f(ff) {}

	double operator()(double x){
		if(f=="tanh0") {
			if(x<=0){
				return v+v*tanh(x/v);
			} else {
				return v+(2-v)*tanh(x/(2-v));
			}
		} else if(f=="tanh1"){
			if(x<=0){
				return v*tanh(x/v);
			} else {
				return (2-v)*tanh(x/(2-v));
			}
		}
		else if(f=="th_lin") {
			if(x<=0) {
				return 0;
			} else {
				return pow(x,v);
			}
		} else {
			if(x<=0){
				return v*tanh(x/v);
			} else {
				return (2-v)*tanh(x/(2-v));
			}
		}

	}
};





#endif
