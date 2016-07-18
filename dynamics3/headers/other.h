#ifndef GUARD_other_h
#define GUARD_other_h

inline int pow2(int x)
{
	if (x < 0) return 0;
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x+1;
}

template<class V>
double mean(const V& vec, int size) 
{
	double ret=0.;
	for(int i=0;i<size;++i) {
		ret += vec[i];
	}
	return ret/size;
}

template<class V>
double var(const V& vec, double m, int size)	
{
	double ret=0.;
	for(int i=0;i<size;++i)
		ret += (vec[i] - m)*(vec[i]-m);
	return ret/size;
}


template<class V>
double var(const V& vec, int size)
{
	double m = mean(vec,size);
	return var(vec,m,size);
}


template<class V>
void add_to_vec( V& vec, double a, int size)
{
	for(int i=0;i<size;++i)
		vec[i] += a;
}

template<class V>
void devide_vec( V& vec, double a, int size)
{
	for(int i=0;i<size;++i) vec[i] /= a;
}

template<class V>
void normal_form(V& vec, double m, double std, int size)
{
	for(int i=0;i<size;++i)
		vec[i] = (vec[i]-m)/std;
}

template<class V>
void normal_form(V& vec, int size)
{
	double m = mean(vec,size);
	double std = sqrt(var(vec,m,size));
	for(int i=0;i<size;++i)
		vec[i] = (vec[i]-m)/std;
}
#endif

