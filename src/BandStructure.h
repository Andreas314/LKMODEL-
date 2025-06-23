#pragma once
#include <vector>
#include <string>
#include <array>
#include "Hamiltonian_functions_real.h"
class band_structure{
public:
	band_structure(hamiltonian &input);
	void do_plot();
private:
	void read_kpoints();
	double string_to_double(std::string val);
	int string_to_int(std::string val);
	void find_begin(std::ifstream &input);
	void fill_kpoints(std::ifstream &input);
	void make_actual_kpath();
	std::array<double, 3> find_difference(std::array<double, 3> &begin, std::array<double, 3> &end);
	std::array<double, 3> make_step(std::array<double, 3> &begin, std::array<double, 3> &step_direction, int num_steps);
	double norm_of(std::array<double, 3> &vec);
	void compute_eigenvals();
	std::array<std::vector<double>, 8> bands;
	std::vector<std::array<double,3>> kpath, kpath_points;;
	std::vector<double> kpath_x_axis;
	int len_kpath, num_kpoints;
	std::vector<std::string> labels;
	hamiltonian MyHamiltonian;

};
