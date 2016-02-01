#ifndef GUARD_functions_h
#define GUARD_functions_h

// transfer functions

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
