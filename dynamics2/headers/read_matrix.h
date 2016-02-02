#ifndef GUARD_read_matrix_h
#define GUARD_read_matrix_h

#include <string>
#include <iostream>
#include <fstream>

template<class V>
void read_matrix(V& vec,int N, std::string file_name,char delimiter = ';')
{
	std::ifstream in;
	try { 
		 in.open(file_name);
	} catch(exception& e) {
		std::cerr << " ERROR in opening file "
		<< file_name << '\t' << e.what() <<'\t';
		return;
	}
	string temp;
	std::getline(in,temp,delimiter);
	std::string::const_iterator it1 = temp.begin();
	std::string::const_iterator it2 = temp.begin();
	std::string::const_iterator b = temp.begin();
	std::string::const_iterator e = temp.end();
	for(int i=0;i<N;++i) {
		try {
			while(*it2!=delimiter and it2!=e) ++it2;
			string dd(temp,it1-b,it2-it1);
			vec[i] = std::stod(dd);
			it1=++it2;
		} catch(exception& e) {
			std::cerr << " ERROR in reading file " << file_name
			<< '\t' << e.what() <<'\n';
			return;
		}
	}
}

template<class M>
void read_matrix(M& mat,int Nrow, int Ncol,std::string file_name,
		char delimiter1 = '\n', char delimiter2 = ';')
{
	std::ifstream in;
	try { 
		 in.open(file_name);
	} catch(exception& e) {
		std::cerr << " ERROR in opening file "
		<< file_name << '\t' << e.what() <<'\t';
		return;
	}
	for(int i=0;i<Nrow;++i) {
		string temp;
		std::getline(in,temp,delimiter1);
		std::string::const_iterator it1 = temp.begin();
		std::string::const_iterator it2 = temp.begin();
		std::string::const_iterator b = temp.begin();
		std::string::const_iterator e = temp.end();
		for(int j=0;j<Ncol;++j) {
			try {
				while(*it2!=delimiter2 and it2!=e) ++it2;
				string dd(temp,it1-b,it2-it1);
				mat[i][j] = std::stod(dd);
			} catch(exception& e) {
				std::cerr << " ERROR in reading file " << file_name
				<< '\t' << e.what() <<'\n';
				return;
			}
		}
	}
}

#endif
