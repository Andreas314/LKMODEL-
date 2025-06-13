#define f_2 1/sqrt(2)
#define f_3 1/sqrt(3)
#define f_6 1/sqrt(6)
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

optical_matrix::optical_matrix(hamiltonian &input): MyHamiltonian(move(input)), px(make_unique<double[]>(4 * 8 * 8)), 
			py(make_unique<double[]>(4 * 8 * 8)), pz(make_unique<double[]>(4 * 8 * 8)),
			px_complex(make_unique<complex<double>[]>(8 * 8)),
			py_complex(make_unique<complex<double>[]>(8 * 8)),
			pz_complex(make_unique<complex<double>[]>(8 * 8)){
	size = MyHamiltonian.get_size();
	assemble_px();
	assemble_py();
	assemble_pz();
}
optical_matrix::~optical_matrix() = default;

void optical_matrix::get_p(int n, unique_ptr<double[]> &input){
	if (n < 0 && n > 3){
		throw runtime_error("Error: Bad index for optical matrices, n has to be between 0 and 2");
	}
	vector<unique_ptr<double[]> *> matrices = {&px, &py, &pz};
	for (int ii = 0; ii < 256; ++ii){
		input[ii] = (*(matrices[n]))[ii];
	}
	if (n == 0){
		assemble_px();
	}
	else if (n == 1){
		assemble_py();
	}
	else if (n == 2){
		assemble_pz();
	}
}

void optical_matrix::get_p_complex(int n, unique_ptr<complex<double>[]> &input){
	if (n < 0 && n > 3){
		throw runtime_error("Error: Bad index for optical matrices, n has to be between 0 and 2");
	}
	vector<unique_ptr<complex<double>[]> *> matrices = {&px_complex, &py_complex, &pz_complex};
	for (int ii = 0; ii< 64; ii++){
		input[ii] = (*(matrices[n]))[ii];
	}
	if (n == 0){
		assemble_px();
	}
	else if (n == 1){
		assemble_py();
	}
	else if (n == 2){
		assemble_pz();
	}
}


void optical_matrix::compute_at_kpoint(double kx, double ky, double kz){
	MyHamiltonian.diagonal_at_k_point(kx, ky, kz);
	compute_p();
}

double optical_matrix::get_energy(int n){
	double energy = MyHamiltonian.get_value(n);
	return energy;
}

void optical_matrix::to_complex_trans(unique_ptr<double[]> &inp, unique_ptr<complex<double>[]> &output){
	int kk = 0, n= 0;
	for (int ii = 0; ii < 16; ++ii){
		if (kk < 2){
			for (int jj = 0; jj < 8; jj++){
				double real_part = inp[ii * 16 + jj];
				double imaginary_part = inp[ii * 16 + 8 + jj];
				output[n * 8 + jj] = real_part + 1i * imaginary_part;
			}
			n++;
		}
		if (kk == 3){
			kk = 0;
			continue;
		}
		kk++;
		}
	}
void optical_matrix::to_complex(unique_ptr<double[]> &inp, unique_ptr<complex<double>[]> &output){
	for (int ii = 0; ii < 8; ++ii){
		for (int jj = 0; jj < 8; jj++){
			double real_part = inp[ii * 16 + jj];
			double imaginary_part = inp[(ii + 8) * 16 + jj];
			output[ii * 8 + jj] = real_part + 1i * imaginary_part;
		}
	}
}
void optical_matrix::compute_p(){
		vector<unique_ptr<double[]> *> matrices = {&px, &py, &pz};
		vector<unique_ptr<complex<double>[]> *> matrices_comp = {&px_complex, &py_complex, &pz_complex};
		unique_ptr<double[]> trans_matrix = make_unique<double[]>(4 * 8 * 8);
		unique_ptr<complex<double>[]> buffer =  make_unique<complex<double>[]>( 8 * 8);
		unique_ptr<complex<double>[]> comp_trans_matrix = make_unique<complex<double>[]>(8 * 8);
		for (auto matrix : matrices_comp){
			MyHamiltonian.make_trans_matrix(trans_matrix);
			to_complex_trans(trans_matrix, comp_trans_matrix);
			Matmul_Wrapper('N', 'N', comp_trans_matrix, *matrix, buffer, 8, 8, 8, 1.0 + 0i);
			Matmul_Wrapper('N', 'C', buffer, comp_trans_matrix, *matrix, 8, 8, 8, 1.0 + 0i);
		}
}

void optical_matrix::copy_buffer(unique_ptr<complex<double>[]> &to, unique_ptr<complex<double>[]> &from){
	for (int ii = 0; ii < 8 * 8; ++ii){
		to[ii] = from[ii];
	}
}

constexpr void optical_matrix::assemble_px(){
	set<int> indices = {9, 14, 15, 24, 44, 60, 74, 75, 77, 92, 104, 120};
	vector<double> values = {-f_2, - f_6, - f_3, f_2, - f_6, - f_3, f_6, f_3, f_2, - f_2, f_6, f_3};
	int counter = 0;
	for (int ii = 0; ii < 2 * 8 * 8; ++ii){
		if (indices.contains(ii)){
			px[ii] = values[counter];
			++counter;
		}
		else{
			px[ii] = 0;
		}
	}
	copy_other_half(px);
	to_complex(px, px_complex);
}

constexpr void optical_matrix::assemble_py(){
	set<int> indices = {1, 6, 7, 16, 36, 52, 66, 67, 69, 84, 96, 112};
	vector<double> values = {-f_2, f_6, f_3, - f_2, f_6, f_3, f_6, f_3, f_2, f_2, f_6, f_3};
	int counter = 0;
	for (int ii = 0; ii < 2 * 8 * 8; ++ii){
		if (indices.contains(ii)){
			py[ii] = values[counter];
			++counter;
		}
		else{
			py[ii] = 0;
		}
	}
	copy_other_half(py);
	to_complex(py, py_complex);
}

constexpr void optical_matrix::assemble_pz(){
	set<int> indices = {10, 11, 40, 56, 78, 79, 108, 124};
	vector<double> values = {2 * f_6, - f_3, - 2 * f_6, f_3, 2 * f_6, - f_3, - 2 * f_6, f_3};
	int counter = 0;
	for (int ii = 0; ii < 2 * 8 * 8; ++ii){
		if (indices.contains(ii)){
			pz[ii] = values[counter];
			++counter;
		}
		else{
			pz[ii] = 0;
		}
	}
	copy_other_half(pz);
	to_complex(pz, pz_complex);
}

void optical_matrix::copy_other_half(unique_ptr<double[]> &values){
	for (int ii = 128; ii < 256; ++ii){
		if ((ii / 8) % 2 == 0){
			values[ii] = -values[ii - 128 + size];
		}
		else{
			values[ii] = values[ii - 128 - size];
		}
	}
}

