#ifndef GUARD_acorr_analysis_h
#define GUARD_acorr_analysis_h

#include "acorr.h"

#include <vector>

void acorr_analysis( const std::vector<std::vector<double> >& xt,
			const std::vector<std::vector<double> >& phixt,
			const std::vector<char> EI,
			VecDoub& deltaE,VecDoub& deltaI,VecDoub& psd_deltaE,
			VecDoub& psd_deltaI, VecDoub& CE, VecDoub& CI, VecDoub& psd_CE,
			VecDoub& psd_CI,int N, int Ne,int Ni,int tsave,int t2)
{
	VecDoub temp_dx(t2,0.0);
	VecDoub temp_phi(t2,0.0);
	VecDoub temp_delta(t2,0.0);
	VecDoub temp_psd_delta(t2,0.0);
	VecDoub temp_C(t2,0.0);
	VecDoub temp_psd_C(t2,0.0);


	//population means
	double pmE = 0.0;
	double pmI = 0.0;
	double pmi;
	for(int i=0;i<N;++i) {
		for(int ti=0;ti<tsave;++ti) {
			if(EI[i]=='E') pmE += xt[i][ti]/(Ne*tsave);
			else pmI += xt[i][ti]/(Ni*tsave);
		}
	}

	for(int i=0;i<N;++i) {
		if(EI[i]=='E') pmi = pmE;
		else pmi = pmI;

		for(int ti=0;ti<tsave;++ti) {
			temp_dx[ti] = xt[i][ti]-pmi;
			temp_phi[ti] = phixt[i][ti];
		}

		if(EI[i]=='E'){
			autocorrel_psd(temp_dx,tsave,temp_delta,temp_psd_delta);
			autocorrel_psd(temp_phi,tsave,temp_C,temp_psd_C);
		} else {
			autocorrel_psd(temp_dx,tsave,temp_delta,temp_psd_delta);
			autocorrel_psd(temp_phi,tsave,temp_C,temp_psd_C);
		}

		for(int ti=0;ti<tsave;++ti) {
			if(EI[i]=='E') {
				deltaE[ti] += temp_delta[ti]/(double)Ne;
				CE[ti] += temp_C[ti]/(double)Ne;
			} else {
				deltaI[ti] += temp_delta[ti]/(double)Ni;
				CI[ti] += temp_C[ti]/(double)Ni;
			}
		}			

		for(int ti=0;ti<t2/2;++ti) {
			if(EI[i]=='E') {
				psd_deltaE[ti] += temp_psd_delta[ti]/(double)Ne;
				psd_CE[ti] += temp_psd_C[ti]/(double)Ne;
			} else {
				psd_deltaI[ti] += temp_psd_delta[ti]/(double)Ni;
				psd_CI[ti] += temp_psd_C[ti]/(double)Ni;
			}
		}	

	}
}


void acorr_analysis(const std::vector<std::vector<double> >& xt,
		const std::vector<std::vector<double> >& phixt,
		VecDoub& delta, VecDoub& psd_delta,VecDoub& C,
		VecDoub& psd_C, int N,int tsave,int t2)
{
	VecDoub temp_dx(t2,0.0);
	VecDoub temp_phi(t2,0.0);
	VecDoub temp_delta(t2,0.0);
	VecDoub temp_psd_delta(t2,0.0);
	VecDoub temp_C(t2,0.0);
	VecDoub temp_psd_C(t2,0.0);


	//population means
	double pm = 0.0;

	for(int i=0;i<N;++i) {
		for(int ti=0;ti<tsave;++ti) {
			pm += xt[i][ti]/(N*tsave);
		}
	}

	for(int i=0;i<N;++i) {
		for(int ti=0;ti<tsave;++ti) {
			temp_dx[ti] = xt[i][ti]-pm;
			temp_phi[ti] = phixt[i][ti];
		}

		autocorrel_psd(temp_dx,tsave,temp_delta,temp_psd_delta);
		autocorrel_psd(temp_phi,tsave,temp_C,temp_psd_C);

		for(int ti=0;ti<tsave;++ti) {
			delta[ti] += temp_delta[ti]/(double)N;
			C[ti] += temp_C[ti]/(double)N;
		}			

		for(int ti=0;ti<t2/2;++ti) {
			psd_delta[ti] += temp_psd_delta[ti]/(double)N;
			psd_C[ti] += temp_psd_C[ti]/(double)N;
		}	

	}

}
#endif