#pragma once

#include <complex>
#include <iostream>
#include <string>
#include <set>
#include <memory>
#include <vector>
#include <map>
using c_d = std::complex<double>;

class hamiltonian{
public:
	hamiltonian();
	hamiltonian(std::string filename_1);
	~hamiltonian();
	hamiltonian(hamiltonian &&other);
	void diagonal_at_k_point(double kx, double ky, double kz);
	void make_trans_matrix(std::unique_ptr<double[]>& array);
	void make_trans_matrix_T(std::unique_ptr<double[]>& array);
	void get_vector(int n, std::unique_ptr<double[]>& array);
       	double get_value(int n);
	std::string get_filename();
	double get_lattice_const();
	int get_size();
	double get_ep();
	void print_ham();
	static std::map<std::string, double> parameters;
private:
	static std::set<std::string> expected_parameters;
	static std::set<std::string> extra_parameters;
	void assemble(double kx, double ky, double kz);
	void to_matrix(c_d (&input)[9]);
	void copy_other_half();
	void copy_to_bottom();
	void read_values(std::string filename_1);
	void tolower_string(std::string &input);
	double string_to_double(std::string val, std::string keyword);
	void make_tracer(std::map<std::string, bool> &input);
	void is_all_loaded(std::map<std::string, bool> &input);
	int size;
	static std::string filename; 
	std::unique_ptr<double[]> values;
	std::unique_ptr<double[]> eigen_vals;
	bool info;

};
