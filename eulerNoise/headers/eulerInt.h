#ifndef GUARD_eulerInt_h
#define GUARD_eulerInt_h

#include "../Eigen/Core"
#include <random>
#include <math.h>
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


//
// euler integration with first order approximation of noise 
// with std of stdN
// each integration step a random vector is added with 
// zero mean and std of sqrt(dt)*stdN
void eulerInt(Eigen::VectorXd& x, Eigen::VectorXd& r,
		const Eigen::MatrixXd& J, const FI& phi, 
		const double& dtsave, const double& dt,
		double stdN, std::default_random_engine& generator,
	    std::normal_distribution<double>& Ndist)
{
	stdN *= sqrt(dt);
	double t =0;
	// integrate until t+dt exceeds dtsave 
	while((t+dt)<=dtsave) {
		x = (1-dt)*x +dt*J*r;
		for(int i=0;i<x.size();++i)
			x(i) += stdN*Ndist(generator);
		phi(r,x);
		t += dt;
	}
	// integrate remaining part
	stdN *= sqrt(dtsave-t)/sqrt(dt);
	if(t<dtsave){
		x = (1-dtsave+t)*x + (dtsave-t)*J*r;
		for(int i=0;i<x.size();++i)
			x(i) += stdN*Ndist(generator);
		phi(r,x);
	}
}



#endif

