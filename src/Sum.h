#include <array>
#include <vector>
#include <complex>
typedef std::array<std::array<std::array<std::complex<double>, 6>,3>,3> tensor;
void write_to_file(std::ofstream &file,auto &last_tensor, double omega);
void sum_and_save(auto &result, auto & current, int num_thr);
void sum_and_save_decomposed(auto &result, auto & current, int num_thr);
std::vector<tensor> make_qi_tensor(int do_decomposition);
void read_omegas(std::array<double, 3>& inp);
double make_sigma(double d_k , double lattice);
double gauss(double x, double x0, double sigma);
