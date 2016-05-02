#ifndef GUARD_eulerInt_h
#define GUARD_eulerInt_h

#include "../Eigen/Core"
#include "functions.h"


void eulerInt(Eigen::VectorXd& x, Eigen::VectorXd& r,
		const Eigen::MatrixXd& J, const FI& phi, 
		const double& dtsave, const double& dt)
{
	double t =0;
	// integrate until t+dt exceeds dtsave 
	while((t+dt)<=dtsave) {
		x = (1-dt)*x +dt*J*r;
		phi(r,x);
		t += dt;
	}
	// integrate remaining part
	if(t<dtsave){
		x = (1-dtsave+t)*x + (dtsave-t)*J*r;
		phi(r,x);
	}
}

#endif

