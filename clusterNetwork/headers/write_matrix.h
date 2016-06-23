#ifndef GUARD_write_matrix_h
#define GUARD_write_matrix_h


#include "../Eigen/Core"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

template<typename Vec>
void write_matrix(const Vec &v, const int size, std::string name,int prec = 16, char delimiter='\n')
{
	std::ofstream outstream(name);
	outstream << std::setprecision(prec);	
	for(int i=0;i<size;i++) outstream << v[i] << delimiter;
}

template<typename Mat>
void write_matrix(const Mat &m, const int cols, const int rows, std::string name,int prec=16,
				char delimiter1=';', char delimiter2='\n')
{
	std::ofstream outstream(name);
	outstream << std::setprecision(prec);
	for(int r=0;r<cols;r++){
		for(int c=0;c<rows;c++) {
			outstream << m[r][c];
			if( c<(rows-1) ) outstream << delimiter1;
		}
		outstream << delimiter2;
	}
}

template<>
void write_matrix<Eigen::MatrixXd>(const Eigen::MatrixXd& m, const int cols,
	const int rows, std::string name, int prec, char delimiter1, char delimiter2)
{
	std::ofstream outstream(name);
	outstream << std::setprecision(prec);
	for(int r=0;r<cols;++r) {
		for(int c=0;c<rows;++c) {
			outstream << m(r,c);
			if(c<(rows-1) ) outstream << delimiter1;
		}
		outstream << delimiter2;
	}
}

#endif

