#ifndef GUARD_read_input_h
#define GUARD_read_input_h

#include <fstream>
#include <string>


std::string remove_text(std::string str)
{
	std::size_t pos = str.find('=');
	return str.substr(pos+1);
}

void read_input(int& N, double& g,int& tf, int& tsave,int& seed,std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	g = std::stod(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	tf = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tsave = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	seed = std::stoi(temp);


}



#endif


