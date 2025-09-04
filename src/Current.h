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

typedef std::array<std::complex<double>, 3> Jvector;
typedef std::array<std::array<std::array<std::complex<double>, 6>,3>,3> tensor;
class Current{
public:
	Current(std::string filename);
	void write_J();
	Jvector to_solve(double t, Jvector y);
private:
	std::map<std::string, double> params;
	std::map<std::string, std::array<double, 3>> params_vector;
	std::vector<tensor> qi_tensor;
	std::vector<double> omegas;
	Jvector factor;
	std::vector<Jvector> J;
	double gauss_t(double t);
	std::complex<double> gauss_omega(double omega);
	void make_factor();
	void solve();
};
