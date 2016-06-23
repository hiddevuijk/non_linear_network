#ifndef GUARD_getStat_h
#define GUARD_getStat_h


void getStat(const Eigen::MatrixXd& xE,
			const Eigen::MatrixXd& xI,
			Eigen::MatrixXd& m0,
			Eigen::MatrixXd& m1,
			Eigen::MatrixXd& m2,
			const int& Nstati)
{
	double temp;
	for(int ri=0;ri<xE.rows();++ri) {
		m0(ri,Nstati) = xE.row(ri).mean();
		m1(ri,Nstati) = 0.0;
		m2(ri,Nstati) = 0.0;
		for( int ci=0;ci<xE.cols();++ci) {
			temp = (xE(ri,ci) - m0(ri,Nstati))*(xE(ri,ci) - m0(ri,Nstati))/xE.cols();
			m1(ri,Nstati) += temp; 
			m2(ri,Nstati) += temp*(xE(ri,ci) - m0(ri,Nstati));
		}
	}
	for(int ri=0;ri<xI.rows();++ri) {
		m0(ri+xE.rows(),Nstati) = xI.row(ri).mean();
		m1(ri+xE.rows(),Nstati) = 0.0;
		m2(ri+xE.rows(),Nstati) = 0.0;
		for( int ci=0;ci<xI.cols();++ci) {
			temp = (xI(ri,ci) - m0(ri+xE.rows(),Nstati))*(xI(ri,ci) - m0(ri+xE.rows(),Nstati))/xI.cols();
			m1(ri+xE.rows(),Nstati) += temp; 
			m2(ri+xE.rows(),Nstati) += temp*(xI(ri,ci) - m0(ri+xE.rows(),Nstati));
		}
	}
	
}



void getStat(const Eigen::MatrixXd& x,
			Eigen::MatrixXd& m0,
			Eigen::MatrixXd& m1,
			Eigen::MatrixXd& m2,
			const int& Nstati)
{
	double temp;
	for(int ri=0;ri<x.rows();++ri) {
		m0(ri,Nstati) = x.row(ri).mean();
		m1(ri,Nstati) = 0.0;
		m2(ri,Nstati) = 0.0;
		for( int ci=0;ci<x.cols();++ci) {
			temp = (x(ri,ci) - m0(ri,Nstati))*(x(ri,ci) - m0(ri,Nstati))/x.cols();
			m1(ri,Nstati) += temp; 
			m2(ri,Nstati) += temp*(x(ri,ci) - m0(ri,Nstati));
		}
	}
}


#endif
