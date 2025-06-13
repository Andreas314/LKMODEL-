#define pi_my 3.1415926
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cstdlib>
#include "Hamiltonian_functions_real.h"
#include "KMesh.h"
using namespace std;
int string_to_int(string &val){
try{
	int value_d;
	value_d = stoi(val);
	return value_d;
	}
catch(invalid_argument){
	throw runtime_error("Error: Bad value "+val+" for KMesh size");
	}
}
double string_to_double(string &val){
try{
	double value_d;
	value_d = stod(val);
	return value_d;
	}
catch(invalid_argument){
	throw runtime_error("Error: Bad value "+val+" for the size of BZ");
	}
}
void tolower_string(string &input){
	for(char &c : input){
		c = tolower(c);
	}
}
void make_kmesh(hamiltonian &Ham, vector<array<double, 3>> &kpoints){
	string filename = Ham.get_filename(), line;
	ifstream open_file(filename);
	string keyword;
	array<int, 3> kmesh_size;
	array<double, 3> bz_part;
	array<array<double, 3>, 3> recciprocal_vectors;
	while (getline(open_file, line)){
		stringstream ss(line);
		ss >> keyword;
		tolower_string(keyword);
		if (keyword == "kmesh"){
			string Ni;
			for (int ii = 0; ii < 3; ++ii){
				ss >> Ni;
				int N = string_to_int(Ni);
				kmesh_size[ii] = N;
			}
		}
		if (keyword == "bz_part"){
			string Ni;
			for (int ii = 0; ii < 3; ++ii){
				ss >> Ni;
				double N = string_to_double(Ni);
				bz_part[ii] = N;
			}
		}
	}
	kpoints.reserve(kmesh_size[0] * kmesh_size[1] * kmesh_size[2]);
	recciprocal_vectors[0] = {bz_part[0], 0, 0};
	recciprocal_vectors[1] = {0 ,bz_part[1], 0};
	recciprocal_vectors[2] = {0, 0, bz_part[2]};
	array<double, 3> point, offset;
	for (int ii = 0; ii < 3; ++ii){
		offset[ii] = 1/kmesh_size[0] * recciprocal_vectors[0][ii];
		offset[ii] += 1/kmesh_size[1] * recciprocal_vectors[1][ii];
		offset[ii] += 1/kmesh_size[2] * recciprocal_vectors[2][ii];
	}
	for (int ii = -kmesh_size[0] /2 ; ii <= kmesh_size[0] / 2; ++ii){
		for (int jj = -kmesh_size[1] / 2; jj <= kmesh_size[1] / 2; ++jj){
			for (int kk = -kmesh_size[2] / 2; kk <= kmesh_size[2] / 2; ++kk){
				for (int hh = 0; hh < 3; ++hh){
					point[hh] = ii / (double)kmesh_size[0] * recciprocal_vectors[0][hh];
					point[hh] += jj / (double)kmesh_size[1] * recciprocal_vectors[1][hh];
					point[hh] += kk / (double)kmesh_size[2] * recciprocal_vectors[2][hh];
					point[hh] += offset[hh];
				}
				kpoints.emplace_back(point);
			}
		}
	}

	
}
