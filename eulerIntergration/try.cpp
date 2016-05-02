#include "headers/write_matrix.h"

#include <iostream>
#include <math.h>
#include <vector>
#include <string>

using namespace::std;

void xtp1(double x, double ti,double tf,double dt) {
	while(ti<tf) {
		x = 1+dt*cos(6*ti);
		ti += dt;
	}
}
int main()
{
	double dt = 0.01;
	vector<double> xt;
	xt.push_back(0);
	double wnew = 0.0;
	for(int i=1;i<xt.size();++i) {
		xt[i] = xtp1(i,dt);
	}
	write_matrix(xt,xt.size(),"xt.csv");

	return 0;
}
