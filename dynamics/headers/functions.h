#ifndef GUARD_functions_h
#define GUARD_functions_h

// transfer functions

double th_linear(double x)
{
	if(x>0) return x;
	else return 0;
}

double tanh01(double x)
{
	
	if(x<=0) return 0.1*tanh(x/0.1);
	else return 1.9*tanh(x/1.9);
}




#endif
