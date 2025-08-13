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

vector<tensor> make_qi_tensor(int do_decomposition){
	// First block to make and open all necessary files
	int val = system("mkdir -p out");
	val = system("touch out/QI_LKMODEL.dat");
	ofstream file("out/QI_LKMODEL.dat");
	ofstream file1, file2, file3;
	if (do_decomposition){
		val = system("touch out/QI_LKMODEL_DECOMP_SO.dat");
		file1.open("out/QI_LKMODEL_DECOMP_SO.dat");
		val = system("touch out/QI_LKMODEL_DECOMP_HH.dat");
		file2.open("out/QI_LKMODEL_DECOMP_HH.dat");
		val = system("touch out/QI_LKMODEL_DECOMP_LH.dat");
		file3.open("out/QI_LKMODEL_DECOMP_LH.dat");
	}
	//write the headers to the files
	file << setw(20) << left << "omega";
	if (do_decomposition){
		file1 << setw(20) << left << "omega";
		file2 << setw(20) << left << "omega";
		file3 << setw(20) << left << "omega";
	}
	array<string, 3> coord_1 = {"x", "y", "z"};
	array<string, 6> coord_2 = {"xx", "yy", "zz", "yz", "xz", "xy"};
	array<string, 2> coord_3 = {"_real", "_imag"};
	for (int aa = 0; aa < 3; aa++){
		for (int cc = 0; cc < 3; cc++){
			for (int bb= 0; bb <  6; bb++){
				for (int dd = 0; dd < 2; dd++){
					string toprint = coord_1[aa] + coord_1[cc] + coord_2[bb] + coord_3[dd];
					file << setw(20) << left << toprint;
					if (do_decomposition){
						file1 << setw(20) << left << toprint;
						file2 << setw(20) << left << toprint;
						file3 << setw(20) << left << toprint;
				}
				}
			}
		}
	}
	file << endl;
	if (do_decomposition){
		file1 << endl;
		file2 << endl;
		file3 << endl;
	}
	
	//read first and last omega and size of the step in eV from INPUT.txt
	array<double, 3> omegas;
	read_omegas(omegas);
	
	//declaration of variables containing the tensor (or its decomposition) for each frequency
	int num_omegas = floor((omegas[1] - omegas[0]) / omegas[2]) + 1;
	vector<tensor> result;
	result.reserve(num_omegas);
	vector<array<tensor, 3>> decomposed_result;
	decomposed_result.reserve(num_omegas);
	
	double omega = omegas[0];
	int ii = 0; //index to govern the loop over frequencies
	
	array<array<int, 2>, 6> symm_ind = {{ {{1,1}}, {{2,2}}, {{3,3}}, {{3,2}}, {{3,1}}, {{2,1}} }};
	array<int, 6> valence_ind = {0,1,2,3,4,5};
	array<int, 2> conduction_ind = {6,7};
	//open hamiltonian, get necessary parameters
	hamiltonian Ham("INPUT.txt");
	double Ep = Ham.get_ep();
	double lattice = Ham.get_lattice_const();
	//make kmesh
	vector<array<double, 3>> kpoints;
	double bz_part = make_kmesh(Ham, kpoints);
	//declare the initial optical matrices objects	
	optical_matrix Opt(Ham);
	//compute the distance between two kpoints
	double delta = pow(kpoints[0][0] - kpoints[1][0], 2);
        delta += pow(kpoints[0][1] - kpoints[1][1], 2);
        delta += pow(kpoints[0][2] - kpoints[1][2], 2);
        delta = sqrt(delta);


	double sigma = 0.05 / (PLANC); 
	double size = kpoints.size();
	//compote the prefactor for the tensor to get it in s^-2mCV^-3
	double Vk = pow( 2 * PI / lattice * delta, 3);
	double P = sqrt(E_MASS * Ep * E_CHARGE / 2);
	double prefactor = Vk / pow(E_MASS, 4) * 2 / pow(PI, 2) * E_CHARGE / pow(PLANC, 3) * pow(P, 4);
	
	//variables for parallel run
	//computed the tensor on each thread for given set of kpoints and then sum it at the end
	int num_thr;
	vector<tensor> current_tensor;
	vector<array<tensor, 3>> decomposed_current_tensor;
	int reserve = 1;
	
	while (omega < omegas[1]){
		#pragma omp parallel shared(omega, do_decomposition, current_tensor, num_thr, kpoints, decomposed_current_tensor) firstprivate(Opt)
		{
		// declaration of the pointers to p matrices
		unique_ptr<complex<double>[]> px = make_unique<complex<double>[]>(64);
		unique_ptr<complex<double>[]> py = make_unique<complex<double>[]>(64);
		unique_ptr<complex<double>[]> pz = make_unique<complex<double>[]>(64);
		#pragma omp master
		{
		num_thr = omp_get_num_threads();
		//reserve the sizes of the current tensors on given frequencies
		//done only once, we only know the size inside the parallel environment 
		if (reserve == 1)
		{
			current_tensor.reserve(num_thr);
			decomposed_current_tensor.reserve(num_thr);
			reserve = 0;
		}
		}
		//fill the tensors with zeros
		int thr_num = omp_get_thread_num();
		for (int a1 = 0; a1 < 3; a1++){
				for (int a2 = 0; a2 < 3; a2++){
					for (int a3 = 0; a3 < 6; a3++){
					current_tensor[thr_num][a1][a2][a3] = 0;
				}
			}
		}
		if (do_decomposition){
		for (int a0 = 0; a0 < 3; a0++){
			for (int a1 = 0; a1 < 3; a1++){
					for (int a2 = 0; a2 < 3; a2++){
						for (int a3 = 0; a3 < 6; a3++){
							decomposed_current_tensor[thr_num][a0][a1][a2][a3] = 0;
						}
					}
				}
			}	
		}
		#pragma omp for
		for (int kk = 0; kk < kpoints.size() ; kk++  ){
			auto kpoint = kpoints[kk];
			//here, the hamiltonian is diagonalized inside the given kpoints (opt matrices are private)
			//p matrices are the computed from the eigenvectors
			//and loaded to px, py and pz
			Opt.compute_at_kpoint(kpoint[0], kpoint[1], kpoint[2]);
			Opt.get_p_complex(0, px);
			Opt.get_p_complex(1, py);
			Opt.get_p_complex(2, pz);
			//declare a vector of pointers to p matrices for easier manipulation
			vector<unique_ptr<complex<double>[]> *> p = {&px, &py, &pz};
			for (int a1 = 0; a1 < 3; a1++){ //loops over coordinates x,y,z
				for (int a2 = 0; a2 < 3; a2++){
					for (int a3 = 0; a3 < 6; a3++){
						int ind1 = symm_ind[a3][0] - 1;
						int ind2 = symm_ind[a3][1] - 1;
						for (int b1 = 0; b1 < 6; b1 ++){ //loop over valence bands
							double valence = valence_ind[b1];
							double E_val = Opt.get_energy(valence);
							for (int b2 = 0; b2 < 2; b2++){ //loop over conduction bands
								double cond = conduction_ind[b2];
								double E_cond = Opt.get_energy(cond);
								double omega_cv = (E_cond - E_val) / (PLANC);
								double omega_cv_1 = (E_cond + E_val) / (PLANC) / 2;
								//compute the weight of each contribution given by
								//replacing the delta function with gaussian
								double weight = gauss(omega_cv, 2 * omega, sigma);
								//the values which are further than 2 sigma from the frequency
								//of the incident light are then cut off for faster run
								if (weight < gauss( -2 * sigma + 2 * omega , 2 * omega, sigma)){
									continue;
								}
								//run over all intermediate states
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
									if (do_decomposition){
										int index = b1 / 2;
										decomposed_current_tensor[thr_num][index][a1][a2][a3] += frac;
									}
								}
							}
						}
					}
				}
			}
		}
		}
		//sum the contribution from each thread and the append it to the resulting list 
		//and write the given element to the file
		sum_and_save(result, current_tensor, num_thr);
		write_to_file(file, result.back(), omega);
		if (file1.is_open() && file2.is_open() && file3.is_open()){
			sum_and_save_decomposed(decomposed_result, decomposed_current_tensor, num_thr);
			array<tensor, 3> last_tensor = decomposed_result.back();
			write_to_file(file1, last_tensor[0], omega);
			write_to_file(file2, last_tensor[2], omega);
			write_to_file(file3, last_tensor[1], omega);
		}
		//here's a place for some nice output (e.g loading bar)
		cout << omega * (PLANC) << ' ' << (result.back())[0][0][0] << endl;
		ii++;
		omega = ii * omegas[2] + omegas[0];
	}
	return result;
}

void write_to_file(ofstream &file,auto &last_tensor, double omega){
	file << setw(10) << setprecision(5) << left <<  omega * (PLANC);
        for (int a1 = 0; a1 < 3; a1++){
            	for (int a2 = 0; a2 < 3; a2++){
                	for (int a3 = 0; a3 < 6; a3++){
				file << setw(20) << setprecision(5) << left << real(last_tensor[a1][a2][a3]);
				file << setw(20) << setprecision(5) << left << imag(last_tensor[a1][a2][a3]);
    			}
		}
	}
	file << endl;
	
}
void sum_and_save(auto &result, auto & current, int num_thr){
	tensor resulting_current_tensor;
	for (int a1 = 0; a1 < 3; a1++){
            	for (int a2 = 0; a2 < 3; a2++){
                	for (int a3 = 0; a3 < 6; a3++){
				for (int thr = 0; thr < num_thr; thr++){
                    			resulting_current_tensor[a1][a2][a3] += current[thr][a1][a2][a3];
				}
    			}
		}
	}
	result.emplace_back(resulting_current_tensor);
}
void sum_and_save_decomposed(auto &result, auto & current, int num_thr){
	array<tensor, 3> resulting_current_tensor;
	for(int a0 = 0; a0 < 3; a0++){
		for (int a1 = 0; a1 < 3; a1++){
            		for (int a2 = 0; a2 < 3; a2++){
                		for (int a3 = 0; a3 < 6; a3++){
					for (int thr = 0; thr < num_thr; thr++){
                    				resulting_current_tensor[a0][a1][a2][a3] += current[thr][a0][a1][a2][a3];
					}
    				}
			}
		}
	}
	result.emplace_back(resulting_current_tensor);
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
