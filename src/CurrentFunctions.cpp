#include <utility>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <complex>
#include <stdexcept>
#include <functional>
#include "CurrentFunctions.h"
using namespace std;
typedef array<complex<double>, 3> Jvector;
typedef array<array<array<complex<double>, 6>,3>,3> tensor;
set<string> expected_parameters = {"i_1", "i_2", "n_1", "n_2", "tau_e", "tau_h", "tau_1", "tau_2", "t_s", "t_f", "dt", "delta_phi", "omega_0", "delay"};
set<string> expected_parameters_vectors = {"e_1", "e_2", "e_3", "j_0"};
set<string> extra_parameters={"begin", "end", "kmesh", "bz_part", "omega_start", "omega_end", "omega_step",
       	"m_c*", "gamma_1^l", "gamma_2^l", "gamma_3^l", "delta_so", "e_g", "e_v", "a", "e_p"};
void read_tensor(string filename, vector<tensor> &qi_tensor, vector<double> &omegas){
	ifstream tensor_file("out/"+filename);
	if (tensor_file.fail()){
		throw runtime_error("Error: "+filename+" not found!");
	}
	string line;
	string dummy;
	getline(tensor_file, line);
	while(getline(tensor_file, line)){
		stringstream ss(line);
		ss >> dummy;
		omegas.push_back(stod(dummy));
		tensor temp_tensor;
		for(int ii = 0; ii < 3; ii++){
			for(int jj = 0; jj < 3; jj++){
				for(int kk = 0; kk < 6; kk++){
					string dummy_r, dummy_c;
					ss >> dummy_r;
					ss >> dummy_c;
					temp_tensor[ii][jj][kk] = stod(dummy_r) + 1i * stod(dummy_c);
				}
			}
		}
		qi_tensor.push_back(temp_tensor);
	}
}
Jvector rg_solve(double tn, Jvector yn, double dt, function<Jvector (double, Jvector)> f){
	Jvector k1 = f(tn, yn);
	for (int ii = 0; ii < 3; ii++){
		k1[ii] *= dt;
	}
	Jvector prep;
	for (int ii = 0; ii < 3; ii++){
		prep[ii] *= yn[ii] + k1[ii] / 2.0;
	}
	

	Jvector k2 = f(tn + dt / 2.0, prep);
	for (int ii = 0; ii < 3; ii++){
		k2[ii] *= dt;
	}
	for (int ii = 0; ii < 3; ii++){
		prep[ii] *= yn[ii] + k2[ii] / 2.0;
	}
	

	Jvector k3 = f(tn + dt / 2.0, prep);
	for (int ii = 0; ii < 3; ii++){
		k3[ii] *= dt;
	}
	for (int ii = 0; ii < 3; ii++){
		prep[ii] *= yn[ii] + k3[ii];
	}
	

	Jvector k4 = f(tn + dt, prep);
	for (int ii = 0; ii < 3; ii++){
		k4[ii] *= dt;
	}


	Jvector y_next;
	for (int ii = 0; ii < 3; ii++){
		y_next[ii] = yn[ii] + 1 / 6.0 * (k1[ii] + 2.0 * k2[ii] + 2.0 * k3[ii] + k4[ii]);
	}
	return y_next;
}

void read_params(map<string, double> &result, map<string, array<double, 3>> &result_vectors){
	ifstream input("INPUT.txt");
	if (input.fail()){
		throw runtime_error("Error: INPUT.txt not found!");
	}
	string line;
	int par_count = 0;
	int par_vect_count = 0;
	while(getline(input, line)){
		stringstream ss(line);
		string dummy;
		ss >> dummy;
		tolower_string_curr(dummy);
		if (dummy == "begin"){
			for (int ii = 0; ii < 5; ii++){
				getline(input, line);
			}
		}
		else if(extra_parameters.contains(dummy) || dummy == "" || dummy[0] == '%'){
			continue;
		}
		else if(expected_parameters.contains(dummy)){
			string dummy_1, dummy_2;
			ss >> dummy_1;
			ss >> dummy_2;
			result[dummy] = stod(dummy_2);
			par_count++;
		}
		else if(expected_parameters_vectors.contains(dummy)){
			array<double, 3> vector;
			for (int ii = 0; ii < 3; ii++){
				string dummy_1;
				ss >> dummy_1;
				vector[ii] = stod(dummy_1);
			}
			result_vectors[dummy] = vector;
			par_vect_count++;
		}
		else{
			throw runtime_error("Error: Unknown parameter "+dummy+"!");
		}
	}
	if (par_vect_count < 4 || par_count < 14){
		throw runtime_error("Error: Too few parameters!");
	}
}

void normal_vect(array<double, 3> &vect){
	double norm = 0;
	for(int ii = 0; ii < 3; ii++){
		norm += vect[ii] * vect[ii];
	}
	norm = sqrt(norm);
	for(int ii = 0; ii < 3; ii++){
		vect[ii] = vect[ii] / norm;
	}
}

void tolower_string_curr(string &input){
	for(char &c : input){
		c = tolower(c);
	}
}

