#ifndef GUARD_read_user_input_h
#define GUARD_read_user_input_h

#include <string>

	// over ride variable if user input is supplied
void read_user_input(int argc, char* argv[], int& N, int& p,
		double& g,double& self, int& seed, int& tf, int& tsave,
		double& tinit,double& r0,double& input_mean,
		double& input_var, double& meanE, double& meanI,
		double& a,int& db,std::string& name)
{
	for(int i=1;i<(argc-1);i+=2) {
		std::string flag = argv[i] ;
		std::string val;
		val = argv[i+1];
		try {val = argv[i+1];}
		catch(exception& e) {}
		if(flag=="-N"){
			try{
				N = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-g"){
			try{
				g = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-self") {
			try {
				self = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<<argv[i] << " not used." <<endl;
			}
		}
		else if(flag=="-a"){
			try{
				a = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-seed"){
			try{
				seed = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-tf"){
			try{
				tf = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-tsave"){
			try{
				tsave = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-tinit"){
			try{
				tinit = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="r0"){
			try{
				r0 = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-p"){
			try{
				int pp = stoi(val);
				if(pp < 3) p = pp;
				else{
					cerr << " error in " << argv[0]
					<< ", "<< argv[i] << " not used."<<endl;
				}
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-meanE"){
			try{
				meanE = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-meanI"){
			try{
				meanI = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-iv"){
			try{
				input_var = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-im"){
			try{
				input_mean = stod(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}


		else if(flag=="-o"){
			try{
				name = val;
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-db") {
			try{
				db = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", " << argv[i] << " not used,"<< endl;
			}
		}

		else {
			cerr << " error in " << argv[0]
			<< ", " << argv[i] << " not used flag."<< endl;
		}
	}

}

void read_user_input(int argc, char* argv[], int& N,
		int& tf, int& tsave,std::string& name)
{
	for(int i=1;i<(argc-1);i+=2) {
		std::string flag = argv[i] ;
		std::string val;
		val = argv[i+1];
		try {val = argv[i+1];}
		catch(exception& e) {}
		if(flag=="-N"){
			try{
				N = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-tf"){
			try{
				tf = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else if(flag=="-tsave"){
			try{
				tsave = stoi(val);
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}


		else if(flag=="-o"){
			try{
				name = val;
			} catch(exception& e) {
				cerr << " error in " << argv[0]
				<< ", "<< argv[i] << " not used."<<endl;
			}
		}

		else {
			cerr << " error in " << argv[0]
			<< ", " << argv[i] << " not used flag."<< endl;
		}
	}

}
#endif
