#ifndef GUARD_read_input_h
#define GUARD_read_input_h

#include <fstream>
#include <string>
#include "functions.h"

std::string remove_text(std::string str)
{
	std::size_t pos = str.find('=');
	return str.substr(pos+1);
}

void read_input(int& N,double& g,double& gm, int& tf, int& tsave,
		int& tinit,double& v, double& h0,double& var_noise,
		int& seed,std::string& name,std::string infile_name)
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
	gm = std::stod(temp);	


	std::getline(input,temp);
	temp = remove_text(temp);
	tf = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tsave = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tinit = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	v = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	h0 = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	var_noise = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	seed = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	name = temp;

}

void read_input(int& N,int& tf,int& tsave,std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	std::getline(input,temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tf = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tsave = std::stoi(temp);

}




void read_integration_vars(double& atol, double& rtol, double& h1, double& hmin, std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	atol = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	rtol = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	h1 = std::stod(temp);


	std::getline(input,temp);
	temp = remove_text(temp);
	hmin = std::stod(temp);

}


#endif


