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

void read_input(int& N,int& p, double& g,int& tf, int& tsave,
		double& r0, double& mean_noise, double& var_noise,
		double& meanE, double& meanI, double& a,int&db,
		int& seed,std::string& name,std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	p = std::stoi(temp);

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
	r0 = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	mean_noise = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	var_noise = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	meanE = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	meanI = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	a = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	db = std::stoi(temp);

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


