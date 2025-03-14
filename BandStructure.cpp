#define tol pow(10,-6)
#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <complex>
#include <cstdlib>
#include "Hamiltonian_functions.h"
#include "BandStructure.h"
using namespace std;

	band_structure::band_structure(hamiltonian &input):MyHamiltonian(move(input)){
		read_kpoints();
	}
	void band_structure::do_plot(){
		system("touch BAND_LKMODEL.dat");
		ofstream open_file("BAND_LKMODEL.dat");
		compute_eigenvals();
		open_file << "\n";
		int len = len_kpath / 2 * num_kpoints;
		int num_bands = MyHamiltonian.get_size();
		for (int ii = 0; ii < num_bands; ++ii){
			for (int jj = 0; jj < len; ++jj){
				open_file << kpath_x_axis[jj] << " ";
				open_file << bands[ii][jj] <<  "\n";
			}
			open_file << "\n";
		}
		
	}
	void band_structure::read_kpoints(){
		string filename_band = MyHamiltonian.get_filename();
		ifstream open_file(filename_band);
		if (open_file.is_open()){
			find_begin(open_file);
			fill_kpoints(open_file);
			make_actual_kpath();
		}
		else{
			throw runtime_error("Error: Input file for band structure not found!");
		}
	}
	double band_structure::string_to_double(string val){
		try{
			double value_d;
			value_d = stod(val);
			return value_d;
		}
		catch(invalid_argument){
			throw runtime_error("Error: Bad value "+val+" for kpath point");
		}
	}
	int band_structure::string_to_int(string val){
		try{
			double value_d;
			value_d = stoi(val);
		return value_d;
		}
		catch(invalid_argument){
			throw runtime_error("Error: Bad value "+val+" for kpath length");
		}
	}
	void band_structure::find_begin(ifstream &input){
		string line, num_points, keyword = "";
		while(keyword != "begin"){
			if(!getline(input, line)){
				throw runtime_error("Error: Kpath not found!");	
			}
			stringstream ss(line);
			ss >> keyword;
			ss >> num_points;
		}
		num_kpoints = string_to_int(num_points);
	}
	void band_structure::fill_kpoints(ifstream &input){
		len_kpath = 0;
		string line;
		bool end_switch = false;
		while(true){
			if(!getline(input, line)){
				throw runtime_error("Error: Wrong kpath format, could reach the end keyword!");
			}
			stringstream ss(line);
			array<double, 3> point;
			for ( int ii = 0; ii < 3; ++ii){
				string value;
				ss >> value;
				if (value == "end"){
					end_switch = true;
					break;
				}
				double coord = string_to_double(value);
				point[ii] = coord;
			}
			if (end_switch){
				break;
			}
			len_kpath++;
			kpath.emplace_back(point);
		}
		if (len_kpath == 0 || len_kpath % 2 != 0){
			throw runtime_error("Error: Number of points in kpath must be nonzero and even!");
		}
	}
	void band_structure::make_actual_kpath(){
		kpath_points.reserve(len_kpath / 2 * num_kpoints);
		kpath_x_axis.reserve(len_kpath / 2 * num_kpoints);
		for (int ii = 0; ii < len_kpath; ii = ii +2){
			array<double, 3> k_begin = kpath[ii];
			array<double, 3> k_end = kpath[ii + 1];
			array<double, 3> difference = find_difference(k_begin, k_end);
			double norm = norm_of(difference);
			for (int jj = 0; jj < num_kpoints; ++jj){
				array<double, 3> result = make_step(k_begin, difference, jj);
				kpath_points.emplace_back(result);
				if (ii == 0 && jj == 0){
					kpath_x_axis.emplace_back(0);
				}
				else if (ii != 0 && jj == 0){
					double last = kpath_x_axis.back();
					kpath_x_axis.emplace_back(last);
				}
				else{
					double last = kpath_x_axis.back();
					kpath_x_axis.emplace_back(last + norm);
				}

			}

		}
	}
	array<double, 3> band_structure::find_difference(array<double, 3> &begin, array<double, 3> &end){
		array<double, 3> diff;
		for (int ii = 0; ii < 3; ++ii){
			diff[ii] = (end[ii] - begin[ii]) / num_kpoints;
		}
		return diff;
	}
	array<double, 3> band_structure::make_step(array<double, 3> &begin, array<double, 3> &step_direction, int num_steps){
		array<double, 3> step;
		for (int ii = 0; ii < 3; ++ii){
			step[ii] = begin[ii] + step_direction[ii] * num_steps;
		}
		return step;
	}
	double band_structure::norm_of(array<double, 3> &vec){
		double norm;
		for (int ii = 0; ii < 3; ++ii){
			norm += vec[ii] * vec[ii];
		}
		return norm;
	}
	void band_structure::compute_eigenvals(){
		for (int ii = 0; ii < 8; ++ii){
			bands[ii].reserve(len_kpath / 2 * num_kpoints);
		}
		int size = MyHamiltonian.get_size();
		int jj = 0;
		for(auto point : kpath_points){
			double kx = point[0];		
			double ky = point[1];		
			double kz = point[2];
			MyHamiltonian.diagonal_at_k_point(kx, ky, kz);
			for (int ii = 0; ii < size; ++ii){
				complex<double> energy = MyHamiltonian.get_value(ii);
				if (imag(energy) > tol){
					cout << "Warning: Imaginary part of the energy in (" << kx <<"," 
						<< ky <<  "," << kz << ") is greater than the tolerance: " << imag(energy) << "!";
				}
				bands[ii].emplace_back(real(energy));	
			}
		}
	}




	/* Variables:
	array<vector<double>, 8> bands
	vector<array<double,3>> kpath, kpath_points
	vector<double> kpath_x_axis;
	int len_kpath, num_kpoints;
	vector<string> labels;
	hamiltonian MyHamiltonian;
	*/

