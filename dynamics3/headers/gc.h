#ifndef GUARD_gc_h
#define GUARD_gc_h

#include <math.h>

double gc(double a, double f)
{
	return 1./sqrt(1-f+f/a);
}


#endif
