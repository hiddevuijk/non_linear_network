#ifndef GUARD_shuffle_h
#define GUARD_shuffle_h

template<class V, class R>
void shuffle(V& vec, int N, R& r)
{
	for(int i=0; i<N;++i) {
		int n = r.int64() % (N -1);
		double temp = vec[i];
		vec[i] = vec[n];
		vec[n] = temp;
	}
}



#endif
