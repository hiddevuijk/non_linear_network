#ifndef GUARD_read_input_h
#define GUARD_read_input_h

#include <fstream>
#include <string>


std::string remove_text(std::string str)
{
	std::size_t pos = str.find('=');
	return str.substr(pos+1);
}

void read_input(int& N, double& meanE,double& meanI,
				double& g,double& a,int& db,int& seed,
				std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	meanE = std::stod(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	meanI = std::stod(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	g = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	a = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	db = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	seed = std::stoi(temp);
}



#endif


