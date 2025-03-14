#define H_PLANC 1.054571817*pow(10, -34)
#define M_ELECTRON 9.1093837*pow(10, -31)
#define J_TO_EV 6.241509*pow(10,18)
#define pi_my 3.1415926
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
#include "Hamiltonian_functions.h"
#include "Wrappers.h"
using namespace std;
using c_d = complex<double>;

map<string,double> hamiltonian::parameters={};
set<string> hamiltonian::expected_parameters = {"m_c*", "gamma_1^l", "gamma_2^l", "gamma_3^l", "delta_so", "e_g", "e_v", "e_c", "a", "e_p"};
set<string> hamiltonian::extra_parameters={"begin", "end"};
string hamiltonian::filename;

hamiltonian::hamiltonian(string filename_1) : size(8),values(make_unique<c_d[]>(size * size)), eigen_vals(make_unique<c_d[]>(size)),eigen_vects(make_unique<c_d[]>(size * size)){
	if (parameters.empty()){
		filename = filename_1;
		read_values(filename_1);
	}
}
hamiltonian::hamiltonian(): size(8),values(make_unique<c_d[]>(size * size)), eigen_vals(make_unique<c_d[]>(size)),eigen_vects(make_unique<c_d[]>(size * size)){
	if (parameters.empty()){
		throw runtime_error("Error: No file to read form!");
	}
}
hamiltonian::~hamiltonian(){}

hamiltonian::hamiltonian(hamiltonian &&other):
	eigen_vals(move(other.eigen_vals)),
	eigen_vects(move(other.eigen_vects)),
	values(move(other.values)){
	size = other.size;}

void hamiltonian::diagonal_at_k_point(double kx, double ky, double kz){
	assemble(kx, ky, kz);
	unique_ptr<c_d[]> dummy = make_unique<c_d[]>(size);
	Eigens_Wrapper(values, eigen_vals, dummy, eigen_vects, size);
}

void hamiltonian::get_vector(int n, unique_ptr<c_d[]>& array){
	for (int ii = 0; ii < 8; ii++){
		array[ii] = eigen_vects[ii * size + n];
	}
}
c_d hamiltonian::get_value(int n){
	return eigen_vals[n];
}

string hamiltonian::get_filename(){
	return filename;
}
double hamiltonian::get_lattice_const(){
	return parameters["a"];
}
int hamiltonian::get_size(){
	return size;
}
void hamiltonian::print_ham(){
	cout << "\n";
	for (int ii = 0; ii < 64; ii++){
		if (ii % 8 == 0 && ii != 0){
			cout << endl;
		}
		printf("|%1.2f + %1.2fi|", real(values[ii]), imag(values[ii]));
		
	}
	cout << "\n";
}

void hamiltonian::assemble(double kx, double ky, double kz){
	for (int ii = 0; ii < size * size; ii++){
		values[ii] = 0;
	}
	// input energies in eV, effective mass dimensionless 
	// and lattice constant in nm
	double E_g = parameters["e_g"];
	double m_eff = parameters["m_c*"];
	double E_v = parameters["e_v"];
	double E_c = parameters["e_c"];
	double a = parameters["a"];
	double delta_so = parameters["delta_so"];
	// input  E_p[eV]
	double E_p = parameters["e_p"];
       	// input k_i dimensionless -> [m^{-1}]
	// therefore a [nm] -> [m]
	kx = (kx * 2 * pi_my / a) * pow(10, 9);
	ky = (ky * 2 * pi_my / a) * pow(10, 9);
	kz = (kz * 2 * pi_my / a) * pow(10, 9);
	//norm [m²]
	double norm_sqr = kx * kx + ky * ky + kz * kz;
	// P0 [sqrt(J² s² eV kg^{-1} (m² m^{-2})] -> P0[sqrt(J eV) m] -> P0[eV m]
	c_d P_0 = -sqrt((H_PLANC * H_PLANC) * E_p / 2 / (M_ELECTRON) * J_TO_EV);
	//dimensionless gammas
	double gamma_1 = parameters["gamma_1^l"] - E_p / (3 * E_g + delta_so);
	double gamma_2 = parameters["gamma_2^l"] - 2 * E_p / (3 * E_g + delta_so);
	double gamma_3 = parameters["gamma_3^l"] - 2 * E_p / (3 * E_g + delta_so);
	double gamma_c = 1/m_eff - E_p / 3 * ( 2 / E_g + 1 / (E_g + delta_so));
	//prefactor in [J² s² kg^{-2}]
	double prefactor = (H_PLANC * H_PLANC) / (2 * M_ELECTRON);
	//Hamiltonian elements in [J]
	double O = prefactor * gamma_c * norm_sqr;
	double P = prefactor * gamma_1 * norm_sqr;
	double Q = prefactor * gamma_2 * (kx * kx + ky * ky - 2 * kz * kz);
	c_d R = prefactor * sqrt(3) * (gamma_2 * (kx * kx - ky * ky) - 2.0 * 1i * gamma_3 * kx * ky);
	c_d S = prefactor * sqrt(6) * gamma_3 * (kx - 1i * ky) * kz;
	c_d T = 1 / sqrt(6) * P_0 * (kx + 1i * kz);
	c_d U = 1 / sqrt(3) * P_0 * kz;
	//Hamiltonian elements in [eV]
	Q *= J_TO_EV;
	P *= J_TO_EV;
	O *= J_TO_EV;
	S *= J_TO_EV;
	R *= J_TO_EV;
	double E_CB = E_c + O;
	double E_HH = E_v - P - Q;
	double E_LH = E_v - P + Q;
	double E_SO = E_v - P - delta_so;
	c_d for_matrix[9] = {E_CB, E_HH, E_LH, E_SO, T, U, S, R, Q};
	to_matrix(for_matrix);
}

void hamiltonian::to_matrix(c_d (&input)[9]){
	for (int ii = 0; ii < 4; ii++){
		values[size * ii + ii] = input[ii];
	}
	values[1] = -sqrt(3) * input[4];
	values[2] = sqrt(2) * input[5];
	values[3] = -input[5];
	values[6] = -conj(input[4]);
	values[7] = -sqrt(2) * conj(input[4]);
	values[8] = -sqrt(3) * conj(input[4]);
	values[10] = sqrt(2) * input[6];
	values[11] = -input[6];
	values[14] = -input[7];
	values[15] = -sqrt(2) * input[7];
	values[16] = sqrt(2) * input[5];
	values[17] = sqrt(2) * conj(input[6]);
	values[19] = -sqrt(2) * input[8];
	values[20] = conj(input[4]);
	values[21] = input[7];
	values[23] = sqrt(3) * input[6];
	values[24] = -input[5];
	values[25] = -conj(input[6]);
	values[26] = -sqrt(2) * input[8];
	values[28] = sqrt(2) * conj(input[4]);
	values[29] = sqrt(2) * input[7];
	values[30] = -sqrt(3) * input[6];
	copy_other_half();
}
void hamiltonian::copy_other_half(){
	for (int ii = 32; ii < 64; ii++){
		if (ii % 8 < 4){
			values[ii+4] = conj(values[ii - 32]);
		}
		else
		{
			values[ii-4] = -conj(values[ii - 32]);
		}
	}
}

void hamiltonian::read_values(string filename_1){
	ifstream open_file(filename_1);
	if (open_file.is_open()){
		string line, keyword, bin, value;
		map<string, bool> tracer;
		make_tracer(tracer);
		while(getline(open_file, line)){
			stringstream ss(line);
			ss >> keyword;
		       	ss >> bin;
			ss >> value;
			tolower_string(keyword);
			if (expected_parameters.contains(keyword)){
				double value_d = string_to_double(value, keyword);
				parameters[keyword] = value_d;
				tracer[keyword] = true;
			}
			else if (keyword == "begin"){
				while (keyword != "end"){
					getline(open_file, line);
					stringstream ss(line);
					ss >> keyword;
				}
			}
			else if (extra_parameters.contains(keyword)){}
			else{
				throw runtime_error("Error: Unknown keyword "+keyword);
				
			}
		}
		is_all_loaded(tracer);
		open_file.close();
	
	}
	else{
		throw runtime_error("Error: Input file not found"); 
	}
}
double hamiltonian::string_to_double(string val, string keyword){
try{
	double value_d;
	value_d = stod(val);
	return value_d;
	}
catch(invalid_argument){
	throw runtime_error("Error: Bad value "+val+" for keyword "+keyword);
	}
}
void hamiltonian::tolower_string(string &input){
	for(char &c : input){
		c = tolower(c);
	}
}

void hamiltonian::make_tracer(map<string, bool> &input){
	for (string key : expected_parameters){
		input[key] = false;
	}
}
void hamiltonian::is_all_loaded(map<string, bool> &input){
	for (auto [key, val] : input){
		if (!input[key]){
			throw runtime_error("Error: Key word missing "+key);
		}
	}
}
