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
#include <iomanip>
#include "Hamiltonian_functions_real.h"
#include "Wrappers.h"
#include "Optical_Matrix.h"
#include "KMesh.h"
#include "Sum.h"
using namespace std;
typedef array<array<array<complex<double>, 6>,3>,3> tensor;

vector<tensor> make_qi_tensor(){
	int val = system("touch QI_LKMODEL.dat");
	ofstream file("QI_LKMODEL.dat");
	file << setw(20) << left << "omega";
	array<string, 3> coord_1 = {"x", "y", "z"};
	array<string, 6> coord_2 = {"xx", "yy", "zz", "yz", "xz", "xy"};
	for (int aa = 0; aa < 3; aa++){
		for (int cc = 0; cc < 3; cc++){
			for (int bb= 0; bb <  6; bb++){
				string toprint = coord_1[aa] + coord_1[cc] + coord_2[bb];
				file << setw(20) << left << toprint;
			}
		}
	}
	file << endl;

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


	double sigma = make_sigma(delta,lattice) / 5;
	double size = kpoints.size();
	//double Vk = 4 * PI / 3 * abs(pow(PI / lattice * bz_part, 3)) / size;
	double Vk = pow(delta * 2 * PI / lattice, 3);
	double P = sqrt(E_MASS * Ep * E_CHARGE / 2);
	double prefactor = Vk / pow(E_MASS, 4) * 2 / pow(PI, 2) * E_CHARGE / pow(PLANC, 3) * pow(P, 4) ;
	complex<double> part_sum = 0;
	int num_thr;
	vector<tensor> current_tensor;
	int reserve = 1;
	while (omega < omegas[1]){
		cout << endl;
		file << setw(20) << setprecision(5) << left <<  omega * (PLANC);
		#pragma omp parallel shared(omega, current_tensor, num_thr, kpoints) firstprivate(Opt, part_sum)
		{
		unique_ptr<complex<double>[]> px = make_unique<complex<double>[]>(64);
		unique_ptr<complex<double>[]> py = make_unique<complex<double>[]>(64);
		unique_ptr<complex<double>[]> pz = make_unique<complex<double>[]>(64);
		#pragma omp master
		{
		num_thr = omp_get_num_threads();
		if (reserve == 1)
		{
			current_tensor.reserve(num_thr);
			reserve = 0;
		}
		}
		int thr_num = omp_get_thread_num();
		for (int a1 = 0; a1 < 3; a1++){ //loop over indicies
				for (int a2 = 0; a2 < 3; a2++){
					for (int a3 = 0; a3 < 6; a3++){
					current_tensor[thr_num][a1][a2][a3] = 0;
				}
			}
		}
		#pragma omp for
		for (int kk = 0; kk < kpoints.size() ; kk++  ){
			auto kpoint = kpoints[kk];
			tensor local_tensor;
			Opt.compute_at_kpoint(kpoint[0], kpoint[1], kpoint[2]);
			Opt.get_p_complex(0, px);
			Opt.get_p_complex(1, py);
			Opt.get_p_complex(2, pz);
			vector<unique_ptr<complex<double>[]> *> p = {&px, &py, &pz};
			for (int a1 = 0; a1 < 3; a1++){ //loop over indicies
				for (int a2 = 0; a2 < 3; a2++){
					for (int a3 = 0; a3 < 6; a3++){
						int ind1 = symm_ind[a3][0] - 1;
						int ind2 = symm_ind[a3][1] - 1;
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
									num *= (*(p[a2]))[cond * 8 + valence];
									complex<double> frac = prefactor * 1i * num / denom * weight;
									current_tensor[thr_num][a1][a2][a3] += frac;
								}
							}
						}
					}
				}
			}
		}
		}
		tensor resulting_tensor;
		for (int a1 = 0; a1 < 3; a1++){
            		for (int a2 = 0; a2 < 3; a2++){
                		for (int a3 = 0; a3 < 6; a3++){
					for (int thr = 0; thr < num_thr; thr++){
                    				resulting_tensor[a1][a2][a3] += current_tensor[thr][a1][a2][a3];
					}
    				}
			}
		}
		result.emplace_back(resulting_tensor);
		cout << endl;
		ii++;
		cout << omega * (PLANC) << ' ' <<  abs(resulting_tensor[0][0][0]) << endl;
        	for (int a1 = 0; a1 < 3; a1++){
            		for (int a2 = 0; a2 < 3; a2++){
                		for (int a3 = 0; a3 < 6; a3++){
					file << setw(21) << setprecision(5) << left << resulting_tensor[a1][a2][a3];
    				}
			}
		}
		file << endl;
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
					inp[0] = stod(value) / (PLANC);
				}
				catch(invalid_argument){
					throw runtime_error("Error: Wrong values for omega_start!");
				}
				counter++;
			}
			if (keyword == "omega_end"){
				try{
					inp[1] = stod(value) / (PLANC);
				}
				catch(invalid_argument){
					throw runtime_error("Error: Wrong values for omega_end!");
				}
				counter++;
			}
			if (keyword == "omega_step"){
				try{
					inp[2] = stod(value) / (PLANC);
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
	return exp( -pow(x - x0, 2) / 2 / pow(sigma, 2) ) /sqrt(2 * PI) / sigma;
	//return 1 / PI * sigma / ( pow(x - x0, 2) + pow(sigma, 2));
}
