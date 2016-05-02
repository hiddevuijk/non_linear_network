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

void read_input(int& N, double& g, double& tf, int& tsave,
		double& tinit,double& dt,int& seed,int& seedI,
		std::string& name,std::string infile_name)
{
	std::ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	g = std::stod(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	tf = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tsave = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tinit = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	dt = std::stod(temp);


	std::getline(input,temp);
	temp = remove_text(temp);
	seed = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	seedI = std::stoi(temp);


	std::getline(input,temp);
	temp = remove_text(temp);
	name = temp;

}

void read_input(int& N,int& tf,int& tsave,std::string infile_name)
{
	std::ifstream input(infile_name);
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
	std::ifstream input(infile_name);
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


