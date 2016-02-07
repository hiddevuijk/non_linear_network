#ifndef GUARD_functions_h
#define GUARD_functions_h

// transfer functions

struct Tanhr0
{
private:
	double r0;
	int p;
public:
	Tanhr0(double r00, int pp) : r0(r00), p(pp) {}

	double operator()(double x){
		if(r0>0) {
			if(p==2) {
				if(x<=0){
					return r0+tanh(x/r0);
				} else {
					return r0+(2-r0)*tanh(x/(2-r0));
				}
			} else {
				if(x<=0){
					return tanh(x/r0);
				} else {
					return (2-r0)*tanh(x/(2-r0));
				}
			}
		} else {
			if(x<=0) return 0;
			else if(x>=1) return 1;
			else return x;
		}
	}

	double derivative(double x) {
		if(r0>0) {
			if(x<=0)
				return 1./(cosh(x/r0)*cosh(x/r0));
			else
				return 1./(cosh(x/(2-r0))*cosh(x/(2-r0)));
		} else {
			if(x < 1 && x >0 )
				return 1.;
			else
				return 0.;
		}
	}
};

double th_linear(double x)
{
	if(x<=0) return 0;
	else if(x>0 and x<1) return x;
	else return 1;
}

double tanh01(double x)
{
	
	if(x<=0) return 0.1*tanh(x/0.1);
	else return 1.9*tanh(x/1.9);
}

double tanh001(double x)
{
	if(x<=0) return 0.1+0.1*tanh(x/0.1);
	else return 0.1+1.9*tanh(x/1.9);
}



#endif
