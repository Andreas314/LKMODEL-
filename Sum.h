#include <array>
#include <vector>
#include <complex>
typedef std::array<std::array<std::array<std::complex<double>, 3>,6>,3> tensor;
std::vector<tensor> make_qi_tensor();
void read_omegas(std::array<double, 3>& inp);
double make_sigma(double d_k , double lattice);
double gauss(double x, double x0, double sigma);
