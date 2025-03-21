#pragma once
#include <memory>
#include "Hamiltonian_functions_real.h"
class optical_matrix{
	public:
		optical_matrix(hamiltonian &Input);
		~optical_matrix();
		void get_px();
	private:
		void compute_p();
		void copy_buffer(std::unique_ptr<double[]> &to, std::unique_ptr<double[]> &from);
		constexpr void assemble_px();
		constexpr void assemble_py();
		constexpr void assemble_pz();
		void copy_other_half(std::unique_ptr<double[]> &values);
		hamiltonian MyHamiltonian;
		int size;
		std::unique_ptr<double[]> px;
		std::unique_ptr<double[]> py;
		std::unique_ptr<double[]> pz;
};
