#ifndef GUARD_functions_h
#define GUARD_functions_h

// transfer functions

struct Tanhr0
{
private:
	double r0;
public:
	Tanhr0(double r00) : r0(r00) {}

	double operator()(double x){
		if(r0>0) {
			if(x<=0){
				return r0+tanh(x/r0);
			} else {
				return r0+(2-r0)*tanh(x/(2-r0));
			}
		} else {
			if(x<=0) return 0;
			else if(x>=1) return 1;
			else return x;
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
