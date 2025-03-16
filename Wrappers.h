#pragma once
#include <memory>
template <typename T>
void Eigens_Wrapper(std::unique_ptr<T[]>& input, std::unique_ptr<T[]>& Eigen_values_real, std::unique_ptr<T[]>& Eigen_values_imaginary, std::unique_ptr<T[]>& Eigen_vectors, int Dimension);
template<typename T>
void Matmul_Wrapper(std::unique_ptr<T[]>& A, std::unique_ptr<T[]>& B, std::unique_ptr<T[]>& result, int Rows_A, int Col_A, int Col_B, T alpha);
template<typename T>
void Eigens_Symm_Wrapper(std::unique_ptr<T[]>& input, std::unique_ptr<T[]>& eigens, int dim);
