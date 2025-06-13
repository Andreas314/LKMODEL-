#define PI 3.1415926
#define E_MASS 9.1093837*pow(10,-31) // kg
#define PLANC 6.5821*pow(10,-16) // eV
#define E_CHARGE 1.60217663*pow(10, -19) // C

#include <memory>
#include <cmath>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <stdexcept>
#include <complex>
#include "Hamiltonian_functions_real.h"
#include "Wrappers.h"
#include "Optical_Matrix.h"
using namespace std;
typedef <array<array<array<double, 3>,6>,3> tensor;

tensor make_qi_tensor(){
	vector<tensor> result;
	array<double, 3> omegas;
	read_omegas(omegas);
	int num_omegas = floor((omega[1] - omega[0]) / omega[2]);
	result.reserve(num_omegas);

	int ii = 0;
	double omega = omegas[0];
	array<array<int, 2>, 6> symm_ind = {{1,1}, {2,2}, {3,3}, {3,2}, {3,1}, {2,1}};
	array<int, 6> valence_ind = {2,3,4,6,7,8};
	array<int, 2> conduction_ind = {1,5};
	
	hamiltonian Ham("INPUT.txt");
	double lat = Ham.get_lattice();
	vector<array<double, 3>> kpoints;
	make_kmesh(Ham, kpoints);	
	optical_matrix Opt(Ham);
	double dk = pow(kpoints[0][0] - kpoints[1][0], 2);
	delta += pow(kpoints[0][1] - kpoints[1][1], 2);
	delta += pow(kpoints[0][2] - kpoints[1][2], 2);
	delta = sqrt(dk);

	unique_ptr<complex<double>[]> px = make_unique<complex<double>[]>(64);
	unique_ptr<complex<double>[]> py = make_unique<complex<double>[]>(64);
	unique_ptr<complex<double>[]> pz = make_unique<complex<double>[]>(64);
	while (omega < omegas[1]){
		tensor current_tensor;
		for (int a1 = 0; a1 < 3; a1++){ //loop over indicies
			for (int a2 = 0; a2 < 6; a2++){
				for (int a3 = 0; a3 < 3; a3++){
					middle_ind_1 = symm_ind[a2][0];
					middle_ind_2 = symm_ind[a2][1];
					tensor[a1][a2][a3] = 0;
					for (auto kpoint : kpoints){
						Opt.compute_at_k_point(kpoint[0], kpoint[1], kpoint[2]);
						Opt.get_p_complex(0, px);
						Opt.get_p_complex(1, py);
						Opt.get_p_complex(2, pz);
						double sigma = make_sigma(delta, kpoint, lattice);
						for (int b1 = 0; b1 < 6; b1 ++){ //loop over bands
							valence = valence_ind[b1];
							E_val = Opt.get_energy[valence];
							for (b2 = 0; b2 < 2; b2++){
								cond = conduction_ind[b2];
								E_cond = Opt.get_energy[cond];
								omega_cv = (E_cond - E_val) / (PLANC);
								omega_cv_1 = (E_cond + E_val) / (PLANC);

								double weight = gauss(omega_cv, omega, sigma);
								for (b3 = 0; b3 < 8; b3++){
									E_inter = Opt.get_energy[b3];
									omega_i = E_inter / (PLANC) * E_CHARGE;
									double denom = omega_cv ^ 3 * (omega_i - omega_cv_1);
									double num = 
								}
							}
						}
					}
				}
			}
		}
		ii++;
		omega += ii * omegas[2];	
	}
		
}

void read_omegas(array<double, 3>& inp){
	ifstream input_file('INPUT.txt');
	if (input_file.is_open()){
		string line;
		int counter;
		while(getline(input_file, line)){
			stringstream ss(line);
			string keyword, bin, value;
			ss >> keyword;
			ss >> bin;
			ss >> value;
			if (keyword == "omega_start"){
				try{
					inp[0] = stod(value);
				}
				catch(invalid_argument){
					throw runtime_error("Error: Wrong values for omega_start!");
				}
				counter++;
			}
			if (keyword == "omega_end"){
				try{
					inp[1] = stod(value);
				}
				catch(invalid_argument){
					throw runtime_error("Error: Wrong values for omega_end!");
				}
				counter++;
			}
			if (keyword == "omega_step"){
				try{
					inp[2] = stod(value);
				}
				catch(invalid_argument){
					throw runtime_error("Error: Wrong values for omega_step!");
				}
				counter++;
			}
		}
		if (counter != 3){
			throw runtime_error("Wrong number of arguemnts for omega range!");
		}
	}
}

double make_sigma(double d_k, double k, double lattice){
	return PLANC / E_MASS * k * 2 * PI * d_k / lattice;
}

double gauss(double x, double x0, double sigma){
	return exp( (x - x0) ^ 2 / 2 / sigma ^ 2 );
}
