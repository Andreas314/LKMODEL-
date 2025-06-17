#define PI 3.1415926
#define E_MASS 9.1093837*pow(10,-31) // kg
#define PLANC 6.5821*pow(10,-16) // eV
#define E_CHARGE 1.60217663*pow(10, -19) // C

#include <memory>
#include <cmath>
#include <utility>
#include <vector>
#include <array>
#include <set>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <complex>
#include <omp.h>
#include "Hamiltonian_functions_real.h"
#include "Wrappers.h"
#include "Optical_Matrix.h"
#include "KMesh.h"
#include "Sum.h"
using namespace std;
typedef array<array<array<complex<double>, 3>,6>,3> tensor;

vector<tensor> make_qi_tensor(){
	vector<tensor> result;
	array<double, 3> omegas;
	read_omegas(omegas);
	int num_omegas = floor((omegas[1] - omegas[0]) / omegas[2]) + 1;
	result.reserve(num_omegas);
	double omega = omegas[0];

	int ii = 0;
	array<array<int, 2>, 6> symm_ind = {{ {{1,1}}, {{2,2}}, {{3,3}}, {{3,2}}, {{3,1}}, {{2,1}} }};
	array<int, 6> valence_ind = {0,1,2,3,4,5};
	array<int, 2> conduction_ind = {6,7};
	
	hamiltonian Ham("INPUT.txt");
	double Ep = Ham.get_ep();
	double lattice = Ham.get_lattice_const();
	vector<array<double, 3>> kpoints;
	double bz_part = make_kmesh(Ham, kpoints);	
	optical_matrix Opt(Ham);
	double delta = pow(kpoints[0][0] - kpoints[1][0], 2);
	delta += pow(kpoints[0][1] - kpoints[1][1], 2);
	delta += pow(kpoints[0][2] - kpoints[1][2], 2);
	delta = sqrt(delta);

	unique_ptr<complex<double>[]> px = make_unique<complex<double>[]>(64);
	unique_ptr<complex<double>[]> py = make_unique<complex<double>[]>(64);
	unique_ptr<complex<double>[]> pz = make_unique<complex<double>[]>(64);

	double sigma = make_sigma(delta,lattice) /5 ;
	double size = kpoints.size();
	double Vk = abs(pow(2 * PI / lattice * bz_part , 3)) / size;
	double P = sqrt(E_MASS * Ep * E_CHARGE / 2);
	double prefactor = Vk / pow(E_MASS, 4) * 2 / pow(PI, 2) * E_CHARGE / pow(PLANC, 3) * pow(P, 4) ;
	complex<double> part_sum = 0;
	while (omega < omegas[1]){
		tensor current_tensor;
		for (int a1 = 0; a1 < 3; a1++){ //loop over indicies
			for (int a2 = 0; a2 < 6; a2++){
				for (int a3 = 0; a3 < 3; a3++){
					current_tensor[a1][a2][a3] = 0;
				}
			}
		}
		#pragma omp parallel for shared(omega, current_tensor) firstprivate(Opt, part_sum)
		for (int kk = 0; kk < kpoints.size(); kk++  ){
			auto kpoint = kpoints[kk];
			Opt.compute_at_kpoint(kpoint[0], kpoint[1], kpoint[2]);
			Opt.get_p_complex(0, px);
			Opt.get_p_complex(1, py);
			Opt.get_p_complex(2, pz);
			vector<unique_ptr<complex<double>[]> *> p = {&px, &py, &pz};
			for (int a1 = 0; a1 < 3; a1++){ //loop over indicies
				for (int a2 = 0; a2 < 6; a2++){
					int ind1 = symm_ind[a2][0] - 1;
					int ind2 = symm_ind[a2][1] - 1;
					for (int a3 = 0; a3 < 3; a3++){
						for (int b1 = 0; b1 < 6; b1 ++){ //loop over bands
							double valence = valence_ind[b1];
							double E_val = Opt.get_energy(valence);
							for (int b2 = 0; b2 < 2; b2++){
								double cond = conduction_ind[b2];
								double E_cond = Opt.get_energy(cond);
								double omega_cv = (E_cond - E_val) / (PLANC);
								double omega_cv_1 = (E_cond + E_val) / (PLANC) / 2;
								double weight = gauss(omega_cv, 2 * omega, sigma);
								for (int b3 = 0; b3 < 8; b3++){
									double E_inter = Opt.get_energy(b3);
									double omega_i = E_inter / (PLANC);
									complex<double> denom = pow(omega_cv , 3) * (omega_i - omega_cv_1);
									complex<double> num = (*(p[a1]))[cond * 8 + cond] - (*(p[a1]))[valence * 8 + valence];
									num *= ((*(p[ind1]))[valence * 8 + b3] * (*(p[ind2]))[b3 * 8 + cond] +
									       	(*(p[ind2]))[valence * 8 + b3] * (*(p[ind1]))[b3 * 8 + cond]) / 2.0;
									num *= (*(p[a3]))[cond * 8 + valence];
									complex<double> frac = prefactor * 1i * num / denom * weight;
									part_sum+=frac;
								}
							}
						}
						current_tensor[a1][a2][a3] = current_tensor[a1][a2][a3] + part_sum;
						part_sum = 0;
					}
				}
			}
		}
		result.emplace_back(current_tensor);
		cout << endl;
		ii++;
		cout << omega * PLANC << ' ' <<  abs(current_tensor[0][0][0]) << endl;
		omega = ii * omegas[2] + omegas[0];
	}
	return result;
}


void read_omegas(array<double, 3>& inp){
	ifstream input_file("INPUT.txt");
	if (input_file.is_open()){
		string line;
		int counter = 0;
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

double make_sigma(double d_k, double lattice){
	return PLANC * E_CHARGE / (E_MASS) * pow(2 * PI, 2) * d_k / pow(lattice,2);
}

double gauss(double x, double x0, double sigma){
	//return exp( -pow(x - x0, 2) / 2 / pow(sigma, 2) ) /sqrt(2 * PI) / sigma;
	return 1 / PI * sigma / ( pow(x - x0, 2) + pow(sigma, 2));
}
