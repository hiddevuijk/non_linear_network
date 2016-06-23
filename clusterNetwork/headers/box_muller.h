#ifndef GUARD_box_muller_h
#define GUARD_box_muller_h

double bm_transform(Ran &r) {
	double x1=0;
	double x2=0;
	double y1=0;
	double y2=0;
	double w=0;
		

	do {
		x1 = 2.0*r.doub() - 1.;
		x2 = 2.0*r.doub() - 1.;
		w = x1*x1 + x2*x2;
	} while( w>= 1.0 || w == 0.0 );

	w = sqrt( (-2.0*log(w))/w);

	return x1*w;
}



#endif

