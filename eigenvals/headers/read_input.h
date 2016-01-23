#ifndef GUARD_read_input_h
#define GUARD_read_input_h

#include <fstream>
#include <string>


std::string remove_text(std::string str)
{
	std::size_t pos = str.find('=');
	return str.substr(pos+1);
}

void read_input(int& N, double& meanS,double& varS, double& meanA, double& varA,int& seed,std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	meanS = std::stod(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	varS = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	meanA = std::stod(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	varA = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	seed = std::stoi(temp);
}



#endif


