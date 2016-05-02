#include "Eigen/Core"

#include <iostream>
#include <vector>
#include <string>

using namespace::Eigen;
using namespace::std;

int main()
{

	VectorXd v(10);
	for(int i=0;i<10;++i) v(i) = 1;
	cout << v << endl;
	VectorXd v1(5);
	v1 = v.head(5);
	cout <<endl << v1 << endl;

	return 0;
}



