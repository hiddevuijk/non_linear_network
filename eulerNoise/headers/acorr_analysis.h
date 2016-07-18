#ifndef GUARD_acorr_analysis_h
#define GUARD_acorr_analysis_h

#include "acorr.h"

#include "../Eigen/Core"
#include <vector>



void acorr_analysis(const Eigen::MatrixXd& xt,
		VecDoub& delta, VecDoub& psd_delta,
		int N,int tsave,int t2)
{
	VecDoub temp_dx(t2,0.0);
	VecDoub temp_delta(t2,0.0);
	VecDoub temp_psd_delta(t2,0.0);


	//population means
	double pm = 0.0;

	for(int i=0;i<N;++i) {
		for(int ti=0;ti<tsave;++ti) {
			pm += xt(i,ti)/(N*tsave);
		}
	}

	for(int i=0;i<N;++i) {
		for(int ti=0;ti<tsave;++ti) {
			temp_dx[ti] = xt(i,ti)-pm;
		}

		autocorrel_psd(temp_dx,tsave,temp_delta,temp_psd_delta);

		for(int ti=0;ti<tsave;++ti) {
			delta[ti] += temp_delta[ti]/(double)N;
		}			

		// normalize with variance
		double varDelta = var_vec(temp_psd_delta,N);
		for(int ti=0;ti<t2/2;++ti) {
			psd_delta[ti] += temp_psd_delta[ti]/(N*varDelta);
		}	

	}

}

#endif
