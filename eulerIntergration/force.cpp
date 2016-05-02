//#include </usr/local/Cellar/eigen/3.2.4/include/eigen3/Eigen/Dense>
#include "Eigen/Core"
// other headers
#include "headers/write_matrix.h"
#include "headers/read_input.h"
#include "headers/read_user_input.h"
#include "headers/step.h"
//std headers
#include <vector>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream>
#include <math.h>
#include <random>

using namespace::std;
using namespace::Eigen;

const double PI = acos(-1);

int main(int argc, char* argv[])
{

	string name;	// output name

	int N;			// number of neurons
	double tInit;	// integrate until tinit befor "real" integration starts
	double tLearnMax;
	double tRun;
	double dt;
	int learnEvery;
	double p;
	double alpha;
	double g;
	int seed;		// seed for the random generator
	
	

	read_input(N,g,tInit,tLearnMax,tRun,dt,learnEvery,p,alpha,seed,name, "input.txt");

	// over ride variable if user input is supplied
	if(argc >1){
		read_user_input(argc,argv,N,g,seed,
			tInit,tLearnMax,tRun,name,dt,p,alpha,learnEvery);
	}

	default_random_engine generator(seed);
	normal_distribution<double> Ndist(0.0,1.0);
	uniform_real_distribution<double> Udist(0.0,1.0);

	if(name!="") name = "_"+name;

	double tMax = tInit+tLearnMax+tRun;
	int tiMax = tMax/dt;
	int tiLearnMax = tLearnMax/dt;
	int tiInitMax =  tInit/dt;
	int tiRunMax = tRun/dt;

	// create connectivity matrix Jgg Jgo W
	MatrixXd Jgg(N,N);
	VectorXd Jgo(N);
	VectorXd W(N);
	VectorXd dW(N);

	// initialize
	for(int i=0;i<N;++i) {
		Jgo(i) = (1-2*Udist(generator));
		W(i) = (1./sqrt(N))*Ndist(generator);
		for(int j=0;j<N;++j) {
			double rNum = Udist(generator);
			if(rNum<p) Jgg(i,j) = g*(1./sqrt(p*N))*Ndist(generator);
		}
	}
	
	VectorXd ft(tiMax);
	VectorXd tt(tiMax);
	double amp = 1.3;
	double freq = 1./60.;
	for(int ti=0;ti<tiMax;++ti){
		ft(ti) = (amp/1.)*sin(1.0*PI*ti*dt*freq);
		ft(ti) += (amp/2.)*sin(2.0*PI*ti*dt*freq);
		ft(ti) += (amp/6.)*sin(3.0*PI*ti*dt*freq);
		ft(ti) += (amp/3.)*sin(4.0*PI*ti*dt*freq);
		ft(ti) /= 1.5;
		tt(ti) = ti*dt;
	}
 

	MatrixXd xt(N,tiMax);
	VectorXd zt(tiMax);
	VectorXd r(N);
	VectorXd err(tiMax);
	VectorXd Wt(tiMax);
	VectorXd WW(tiMax);
	VectorXd k(N);
	MatrixXd P(N,N);


	for(int i=0;i<N;++i) {
		P(i,i) = 1./alpha;
		xt(i,0) = Udist(generator);
	}
	r = xt.col(0).unaryExpr<double(*)(double)>(&tanh);
	zt(0) = W.dot(r);
	err(0) = zt(0) - ft(0);
	Wt(0) = 0.;
	WW(0) = 0.0;

	cout << "Start simulation ... " << endl;
//	double eta = 2.e-3;
	int tiTot=1;
	// start running network, no learning
	for(int ti=1;ti<tiInitMax;++ti){
		xt.col(tiTot) = (1.0-dt)*xt.col(tiTot-1)+Jgg*r*dt + Jgo*zt(tiTot-1)*dt;
		r = xt.col(tiTot).unaryExpr<double(*)(double)>(&tanh);
		zt(tiTot) = r.transpose()*W;
		err(tiTot) = zt(tiTot) - ft(tiTot);
		Wt(tiTot) = Wt(0);
		WW(tiTot) = WW(0);
		++tiTot;
	}

	cout << "Start learning ... " << endl;
	// start learning
	double c;
	for(int ti=1;ti<tiLearnMax;++ti) {
		xt.col(tiTot) = (1.0-dt)*xt.col(tiTot-1)+Jgg*r*dt + Jgo*zt(tiTot-1)*dt;
		r = xt.col(tiTot).unaryExpr<double(*)(double)>(&tanh);
		zt(tiTot) = W.transpose()*r;
		err(tiTot) = zt(tiTot) - ft(tiTot);
		if((ti%learnEvery)==0){
			k = P*r;
			c = 1./(1.+r.transpose()*k);
			P = P - k*k.transpose()*c;
			dW = -1*err(tiTot)*k*c;
			W = W + dW;
			WW(tiTot) = sqrt(W.transpose()*W);
			Wt(tiTot) = sqrt(dW.transpose()*dW);
		} else {
			Wt(tiTot) = Wt(tiTot-1);
			WW(tiTot) = WW(tiTot-1);
		}	
//		eta  = dt*eta*(pow(err(tiTot),1.5) - eta);
//		dW = -1*eta*err(tiTot)*r;
//		Wt(tiTot) = sqrt(dW.dot(dW));
//		W += dW;

		++tiTot;
	}

	cout << "Start post learning simulation ... " << endl;
	// post learning
	for(int ti=1;ti<tiRunMax;++ti) {
		xt.col(tiTot) = (1.0-dt)*xt.col(tiTot-1)+Jgg*r*dt + Jgo*zt(tiTot-1)*dt;
		r = xt.col(tiTot).unaryExpr<double(*)(double)>(&tanh);
		zt(tiTot) = W.transpose()*r;
		err(tiTot) = zt(tiTot) - ft(tiTot);
		Wt(tiTot) = 0;
		WW(tiTot) = 0;
		++tiTot;
	}

	
	write_matrix(zt,zt.size(),"z.csv");
	write_matrix(ft,ft.size(),"f.csv");
	write_matrix(err,err.size(),"e.csv");
	write_matrix(Wt,Wt.size(),"wt.csv");
	write_matrix(WW,WW.size(),"ww.csv");
	return 0;
}


