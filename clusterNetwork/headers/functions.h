#ifndef GUARD_functions_h
#define GUARD_functions_h

#include <string>
#include <math.h>
#include "../Eigen/Core"
// transfer functions

struct FI {
private:
	double r0;
	int p;
public:
	FI(double r00, int pp) : r0(r00), p(pp) {}

	double operator()(double x) const {
		if(p==2) {
			if(x<=0){
				return r0+r0*tanh(x/r0);
			} else {
				return r0+(2-r0)*tanh(x/(2-r0));
			}
		} else {
			if(x<=0){
				return r0*tanh(x/r0);
			} else {
				return (2-r0)*tanh(x/(2-r0));
			}
		}

	}

	void operator()(Eigen::VectorXd& r,const Eigen::VectorXd& x) const {
		if(p==2) {
			for(int i=0;i<x.size();++i) {
				if(x(i)<=0){
					r(i) =  r0+r0*tanh(x(i)/r0);
				} else {
					r(i) =  r0+(2-r0)*tanh(x(i)/(2-r0));
				}
			}
		} else {
			for(int i=0;i<x.size();++i) {
				if(x(i)<=0){
					r(i) = r0*tanh(x(i)/r0);
				} else {
					r(i) = (2-r0)*tanh(x(i)/(2-r0));
				}
			}
		}

	}
};


#endif
