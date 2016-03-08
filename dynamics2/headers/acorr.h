#ifndef GUARD_acorr_h
#define GUARD_acorr_h

/*
	functions for autocorrelation
	both statistical and signal processing
	type. (normalized unnormalized respectively)	

	IMPORTANT!
	depends on nr3.h and fourier.h
	from numerical recipes
*/

#include <math.h>


double mean_vec(const VecDoub& vec, const int N)
{
	double mean=0;
	for(int i=0;i<N;++i) mean += vec[i];
	return mean/N;
}

double var_vec(const VecDoub& vec, const int N, const double& mean)
{
	double var=0;
	for(int i=0;i<N;++i) var += (vec[i]-mean)*(vec[i]-mean);
	return var/N;
}

double var_vec(const VecDoub& vec, const int N)
{
	double m = mean_vec(vec,N);
	return var_vec(vec,N,m);
}

/*
	input: data, zero-padded and size(data)=N=power of 2
	output: ans, size(ans) = N , only non-zero-padded i
			useful. 
	equation:
		C(t_lag) = 1/(n-t_lag) *1/var(data) *E[(data(t)-mean(data))*(data(t+t_lag)-mean(data))]
		With E[] average over t
	example:
			size(data) = n;
			zero-pad dat to size N s.t. N-n is number of time-lags,
			and N a power of 2.

*/
void autocorrel_norm(const VecDoub& data,const int N, VecDoub& ans)
{
	int no2,i,n=data.size();
	double mean = mean_vec(data,N);
	double var = var_vec(data,N,mean);
	for(i=0;i<N;++i) {
		ans[i]=(data[i]-mean)/sqrt(var);
	}
	for(int i=N;i<n;++i) ans[i]=0;
	realft(ans,1);
	no2=n>>1;
	for (i=2;i<n;i+=2) {
		ans[i] = (ans[i]*ans[i]+ans[i+1]*ans[i+1])/no2;
		ans[i+1] = 0;
	}
	ans[0]=ans[0]*ans[0]/no2;
	ans[1]=ans[1]*ans[1]/no2;
	realft(ans,-1);
	for(int i=0;i<N;++i) ans[i]/=(N-i);	
}


/*
	acorr_psd:
		calculates the one sided PSD
		and the autocorrelation(signal processing kind
		so not normalized and not mean subtracted.)
	data:
		Input data. Size(n) must be a power of 2.
		Zero-padded , with # of 0s = number of
		lags.
	N:
		Number of data points, so all point without
		the zero-padded ones.
	ans_acorr:
		The autocorrelation. Size is size of data.
		Only use the first 
	ans_psd:
		The one-sided psd. Size: half of data

*/

void autocorrel_psd(const VecDoub& data, const int N, VecDoub& ans_acorr, VecDoub& ans_psd)
{
	int no2;			// used for normalization of inv. FT
	int i;				// used for loops
	int n=data.size();	// size of data, must be a pow. of 2.

	// copy data
	for(int i=0;i<N;++i) ans_acorr[i] = data[i];
	// set others to 0
	for(int i=N;i<n;++i) ans_acorr[i] = 0;
	// FT data
	realft(ans_acorr,1);

	// Normalization const no2 = n/2
	no2=n>>1;

	// calculation of psd=|X|^2, and
	for (i=2;i<n;i+=2) {
		// psd 
		ans_psd[i/2] =  ans_acorr[i]*ans_acorr[i]+ans_acorr[i+1]*ans_acorr[i+1];
		// /no2 for normalization
		ans_acorr[i] = ans_psd[i/2]/no2;
		// *2 because it is the one-sided psd,
		// /n for normalization s.t. : sum(|x|*2)=sum(psd)
		ans_psd[i/2] *= 2./n;
		// no imag. part in X*compl.conj.(X), because x symmetric
		ans_acorr[i+1]=0;
	}

	// first and last of the psd
	ans_psd[0] = 2*ans_acorr[0]*ans_acorr[0]/n;
	ans_psd[n/2-1] = 2*ans_acorr[1]*ans_acorr[1]/n;
	
	// first and last of the autocorr.
	ans_acorr[0]=ans_acorr[0]*ans_acorr[0]/no2;
	ans_acorr[1]=ans_acorr[1]*ans_acorr[1]/no2;

	// inv. FT to get the autocorrelation
	realft(ans_acorr,-1);

	// devide acorr by number of used points
	for(int i=0;i<N;++i) ans_acorr[i]/=(N-i);
}


/*
	psd_frequencies:
		frequncies of corresponding
		to the psd (one sided)
	N: total number of points in |X|^2
		only N/2 freqs return, because of
		the one sided psd
	dt:
		samplig freq
	freq:
		size N/2
*/
template<class V>
void psd_frequencies(const int N, const double dt, V& freq)
{
	for(int i=0;i<N/2;++i)
		freq[i] = i/(double)(N*dt);
}


/*
	Straightforward implementation
	of the autocorrelation.	

	The signal proccesing kind, so
	not normalized or mean subtracted.
*/

void autocorrel_sf(const VecDoub& data, const int N, VecDoub& ans)
{
	for(int i=0;i<N;i++) {
		for(int j=0;j<(N-i);++j) {
			ans[i] += (data[j+i]*data[j])/(N-i);
		}
	}
}





#endif
