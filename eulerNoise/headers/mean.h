#ifndef GUARD_mean_h
#define GUARD_mean_h

void mean(const Eigen::MatrixXd& xt, VecDoub& avg)
{
	int r = xt.rows();
	int	c = xt.cols();

	for(int ci=0;ci<c;++ci) {
		avg[ci] = 0.;
		for(int ri=0;ri<r;++ri) {
			avg[ci] += xt(ri,ci);
		}
		avg[ci]/=r;
	}
}


#endif
