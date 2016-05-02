#ifndef GUARD_normalize_h
#define GUARD_normalize_h

template<class V>
void normalize(V& v,const int N)
{
	for(int i=N-1;i>=0;--i)
		v[i] /= v[0];
}
		


#endif
