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

void read_tensor(std::string filename, std::vector<tensor> &qi_tensor, std::vector<double> &omegas);
Jvector rg_solve(double tn, Jvector yn, double dt, std::function<Jvector (double, Jvector)> f);
void read_params(std::map<std::string, double> &result, std::map<std::string, std::array<double, 3>> &result_vectors);
void normal_vect(std::array<double, 3> &vect);
void tolower_string_curr(std::string &input);
