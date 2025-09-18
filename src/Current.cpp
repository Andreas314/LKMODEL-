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
#include <functional>
#include <iomanip>
#include "CurrentFunctions.h"
#include "Current.h"
#define pi 3.1415926
#define light_c 3 * pow(10,8)
#define epsilon 8.854 * pow(10, -12)
#define H_PLANC 1.054571817*pow(10, -34)
using namespace std;
typedef array<complex<double>, 3> Jvector;
Current::Current(string filename){
	read_params(params, params_vector);
	read_tensor(filename,qi_tensor, omegas);
}

double Current::gauss_t(double t){
	double tau_1 = params["tau_1"];
	double tau_2 = params["tau_2"];
	double delay = params["delay"];
	double to_return = exp(- t * t / tau_1 / tau_1) * exp(- (t - delay) * (t - delay) / tau_2 / tau_2 / 2);
	return to_return;
}

complex<double> Current::gauss_omega(double omega){
	double tau_1 = params["tau_1"];
	double tau_2 = params["tau_2"];
	double delay = params["delay"];
	complex<double> to_return = exp(- omega * omega * tau_1 * tau_1) * exp(- 4 * omega * omega * tau_2 * tau_2 / 2) * pow(2 * pi, 3 / 2) * tau_2 * tau_1 * tau_1 * exp(- 1i * 2.0 * omega * delay);
	return to_return;
}

void Current::make_factor(){
	array<double, 3> E1 = params_vector["e_1"];
	normal_vect(E1);
	array<double, 3> E2 = params_vector["e_2"];
	normal_vect(E2);
	array<double, 3> E3 = params_vector["e_3"];
	normal_vect(E3);
	double I1 = params["i_1"];
	double I2 = params["i_2"];
	double n1 = params["n_1"];
	double n2 = params["n_2"];
	double omega = params["omega_0"];
	double dpsi = params["delta_psi"];
	int ii = 0;
	double omega_o = omegas[0];
	double omega_n = omegas[1];
	while (omega > omega_o && omega < omega_n){
		try{
			ii++;
			double omega_o = omegas[ii];
			double omega_n = omegas[ii+1];
		}
		catch(const out_of_range oo){
			throw runtime_error("Error: omega_0 not in the computed range!");
		}
	}
	map<int, int> indecies = {{1,1}, {4, 2}, {9, 3}, {6, 4}, {3, 5}, {2, 6}};
	for(int ll = 0; ll < 3; ll++){
		for(int jj = 0; jj < 3; jj++){
			for(int kk = 0; kk < 3; kk++){
				for(int hh = 0; hh < 3; hh++){
					int index = indecies[kk * hh];		
					complex<double> stall = qi_tensor[ii][ll][jj][index] * E1[ll] * E2[ll] * E3[ll] * exp(1i * dpsi);
					factor[ll] += (stall + conj(stall)) * I1 * sqrt(I2) * n1 * sqrt(n2) * pow(1 / (light_c * epsilon), 3.0 / 2.0) * pow(2, 5.0 / 2.0);
				}
			}
		}
	
	}
	cout << endl;
}

Jvector Current::to_solve(double t, Jvector y){
	double delta_phi = params["delta_phi"];
	double tau_h = params["tau_h"];
	double tau_e = params["tau_e"];
	Jvector dydt;
	for (int ii = 0; ii < 3; ii++){
		dydt[ii] = factor[ii] * gauss_t(t)- y[ii] / tau_h * pow(10, 9) - y[ii] / tau_e * pow(10, 9);
	}
	return dydt;
}

void Current::solve(){
	factor = {0, 0, 0};
	make_factor();
	double t1 = params["t_s"];
	double t2 = params["t_f"];
	double dt = params["dt"];
	double len = (t2 - t1) / dt;
	J.reserve(len);
	Jvector J_now;
	for(int ii = 0; ii < 3; ii++){	
		J_now[ii] = params_vector["j_0"][ii];
	}
	double t = t1;
	for(int ii = 0; ii <= len; ii++){
		J.emplace_back(J_now);
		t += dt;
		J_now = rg_solve(t, J_now, dt * pow(10, -9), bind(&Current::to_solve, this, t, J_now));
	}
}

void Current::write_J(){
	solve();
	int dumm = system("touch out/J.dat");
	ofstream output("out/J.dat");
	double t1 = params["t_s"];
	double t2 = params["t_f"];
	double dt = params["dt"];
	double len = (t2 - t1) / dt;
	double t = t1;
	output << setw(20) << left << "t";
	output << setw(20) << left << "J_x";
	output << setw(20) << left << "J_y";
	output << setw(20) << left << "J_z";
	output << "\n";
	for(int ii = 0; ii <= len; ii++){
		output << setw(20) << setprecision(5) << t;
		output << setw(20) << setprecision(5) << real(J[ii][0]);
		output << setw(20) << setprecision(5) << real(J[ii][1]);
		output << setw(20) << setprecision(5) << real(J[ii][2]);
		output << "\n";
		t += dt;
	}
}
