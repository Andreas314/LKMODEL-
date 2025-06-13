#include <iostream>
#include "Wrappers.h"
#include <complex>
#include <memory>
using namespace std;
/*
 * Wrappers around Lapack functions for easier use during calculations
 * Written for 4 types - float, double, complex<float> and complex<double>
 * So far: diagonalization for general matrices
 * TODO: Use lwork in forloops for better performance and optimization of work array
 * TODO: Handle info errors
 * TODO: 
 */

// Eigenvector solver definitions
extern "C" int dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
extern "C" int sgeev_(char*, char*, int*, float*, int*, float*, float*, float*, int*, float*, int*, float*, int*, int*);
extern "C" int cgeev_(char*, char*, int*, complex<float>*, int*, complex<float>*, complex<float>*,  int*, complex<float>*, int*, complex<float>*, int*, float*, int*);
extern "C" int zgeev_(char*, char*, int*, complex<double>*, int*, complex<double>*, complex<double>*,  int*, complex<double>*, int*, complex<double>*, int*, double*, int*);
// Matrix multiplication definition
extern "C" int dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
extern "C" int sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);
extern "C" int cgemm_(char*, char*, int*, int*, int*, complex<float>*, complex<float>*, int*, complex<float>*, int*, complex<float>*, complex<float>*, int*);
extern "C" int zgemm_(char*, char*, int*, int*, int*, complex<double>*, complex<double>*, int*, complex<double>*, int*, complex<double>*, complex<double>*, int*);
// Eigenvector solver for symm matrices
extern "C" int ssyev_(char *, char *, int*, float*, int*, float*, float*, int*, int*);
extern "C" int dsyev_(char *, char *, int*, double*, int*, double*, double*, int*, int*);

//Implementation of Eigen solver wrapper
template<typename T>
void Eigens_Wrapper(unique_ptr<T[]>& input, unique_ptr<T[]>& Eigen_values_real, unique_ptr<T[]>& Eigen_values_imaginary, unique_ptr<T[]>& Eigen_vectors, int Dimension){}
//Implementation for doubles
template <>
void Eigens_Wrapper(unique_ptr<double[]>& input, unique_ptr<double[]>& Eigen_values_real,
	       	unique_ptr<double[]>& Eigen_values_imaginary, unique_ptr<double[]> & Eigen_vectors, int Dimension){
    int info, lwork;
    lwork = 4 * Dimension;//size of the buffer work
    unique_ptr<double[]> work = make_unique<double[]>(lwork); 
    double *EMPTY; //Placeholder for right eigenvectors
    char Eigen_Switch = 'N'; //Right eigenvectors not need
    char Eigen_Switch_On = 'V'; //Left eigenvectors neede
    int dim = 1; // Dimension of right eigenvectors (unused)
    dgeev_(&Eigen_Switch_On, &Eigen_Switch, &Dimension, input.get(), &Dimension, Eigen_values_real.get(),
		    Eigen_values_imaginary.get(), Eigen_vectors.get(), &Dimension, EMPTY, &dim, work.get(), &lwork, &info);
    cout << endl << info << endl;
}
//Implementation for floats
template<>
void Eigens_Wrapper(unique_ptr<float[]>& input, unique_ptr<float[]>& Eigen_values_real,
	       	unique_ptr<float[]>& Eigen_values_imaginary, unique_ptr<float[]> & Eigen_vectors, int Dimension){
    int info, lwork;
    lwork = 4 * Dimension;//size of the buffer work
    unique_ptr<float[]> work = make_unique<float[]>(lwork); 
    float *EMPTY; //Placeholder for right eigenvectors
    char Eigen_Switch = 'N'; //Right eigenvectors not need
    char Eigen_Switch_On = 'V'; //Left eigenvectors neede
    int dim = 1; // Dimension of right eigenvectors (unused)
    sgeev_(&Eigen_Switch_On, &Eigen_Switch, &Dimension, input.get(), &Dimension, Eigen_values_real.get(),
		    Eigen_values_imaginary.get(), Eigen_vectors.get(), &Dimension, EMPTY, &dim, work.get(), &lwork, &info);
    cout << endl << info << endl;
}
//Implementation for complex<float>, BEWARE: Imaginary part of Eigen values is unused, only Eigen_values_real is modified due to use of complex numbers
template<>
void Eigens_Wrapper(unique_ptr<complex<float>[]>& input, unique_ptr<complex<float>[]>& Eigen_values_real,
	       	unique_ptr<complex<float>[]>& Eigen_values_imaginary, unique_ptr<complex<float>[]> & Eigen_vectors, int Dimension){
    int info, lwork;
    lwork = 4 * Dimension;//size of the buffer work
    unique_ptr<complex<float>[]> work = make_unique<complex<float>[]>(lwork); 
    complex<float> *EMPTY; //Placeholder for right eigenvectors
    char Eigen_Switch = 'N'; //Right eigenvectors not need
    char Eigen_Switch_On = 'V'; //Left eigenvectors neede
    int dim = 1; // Dimension of right eigenvectors (unused)
    unique_ptr<float[]> rwork = make_unique<float[]>(2 * Dimension); 
    cgeev_(&Eigen_Switch_On, &Eigen_Switch, &Dimension, input.get(), &Dimension, Eigen_values_real.get(),
		    Eigen_vectors.get(), &Dimension, EMPTY, &dim, work.get(), &lwork, rwork.get(), &info);
    cout << endl << info << endl;
}
//Implementation for complex<double>, BEWARE: Imaginary part of Eigen values is unused, only Eigen_values_real is modified due to use of complex numbers
template<>
void Eigens_Wrapper(unique_ptr<complex<double>[]>& input, unique_ptr<complex<double>[]>& Eigen_values_real,
	       	unique_ptr<complex<double>[]>& Eigen_values_imaginary, unique_ptr<complex<double>[]> & Eigen_vectors, int Dimension){
    int info, lwork;
    lwork = 4 * Dimension; //size of the buffer work
    unique_ptr<complex<double>[]> work = make_unique<complex<double>[]>(lwork); 
    complex<double> *EMPTY; //Placeholder for right eigenvectors
    char Eigen_Switch = 'N'; //Right eigenvectors not need
    char Eigen_Switch_On = 'V'; //Left eigenvectors neede
    int dim = 1; // Dimension of right eigenvectors (unused)
    unique_ptr<double[]> rwork = make_unique<double[]>(2 * Dimension); 
    zgeev_(&Eigen_Switch_On, &Eigen_Switch, &Dimension, input.get(), &Dimension, Eigen_values_real.get(),
		    Eigen_vectors.get(), &Dimension, EMPTY, &dim, work.get(), &lwork, rwork.get(), &info);
}
//Implementation of matrix multiplication wrapper
template<typename T>
void Matmul_Wrapper(char Switch_trans_1, char Switch_trans_2, unique_ptr<T[]>& A, unique_ptr<T[]>& B, unique_ptr<T[]>& result, int Rows_A, int Col_A, int Col_B, T alpha){}
//implementation for double
template<>
void Matmul_Wrapper(char Switch_trans_1, char Switch_trans_2, unique_ptr<double[]>& A, unique_ptr<double[]>& B, unique_ptr<double[]>& result,int Rows_A, int Col_A, int Col_B, double alpha){
	double beta = 0.0;
	dgemm_(&Switch_trans_1, &Switch_trans_2, &Rows_A, &Col_B, &Col_A, &alpha, A.get(), &Rows_A, B.get(), &Col_A, &beta, result.get(), &Rows_A);
}
//implementation for float
template<>
void Matmul_Wrapper(char Switch_trans_1, char Switch_trans_2, unique_ptr<float[]>& A, unique_ptr<float[]>& B, unique_ptr<float[]>& result,int Rows_A, int Col_A, int Col_B, float alpha){
	float beta = 0.0;
	sgemm_(&Switch_trans_1, &Switch_trans_2, &Rows_A, &Col_B, &Col_A, &alpha, A.get(), &Rows_A, B.get(), &Col_A, &beta, result.get(), &Rows_A);
}
//implementation for complex<float>
template<>
void Matmul_Wrapper(char Switch_trans_1, char Switch_trans_2, unique_ptr<complex<float>[]>& A, unique_ptr<complex<float>[]>& B, unique_ptr<complex<float>[]>& result,int Rows_A, int Col_A, int Col_B, complex<float> alpha){
	complex<float> beta = 0.0;
	cgemm_(&Switch_trans_1, &Switch_trans_2, &Rows_A, &Col_B, &Col_A, &alpha, A.get(), &Rows_A, B.get(), &Col_A, &beta, result.get(), &Rows_A);
}
//implementation for complex<double>
template<>
void Matmul_Wrapper(char Switch_trans_1, char Switch_trans_2, unique_ptr<complex<double>[]>& A, unique_ptr<complex<double>[]>& B, unique_ptr<complex<double>[]>& result,int Rows_A, int Col_A, int Col_B, complex<double> alpha){
	complex<double> beta = 0.0;
	zgemm_(&Switch_trans_1, &Switch_trans_2, &Rows_A, &Col_B, &Col_A, &alpha, A.get(), &Rows_A, B.get(), &Col_A, &beta, result.get(), &Rows_A);
}
template<typename T>
void Eigens_Symm_Wrapper(std::unique_ptr<T[]>& input, std::unique_ptr<T[]>& eigens, int dim);

template<>
void Eigens_Symm_Wrapper(std::unique_ptr<float[]>& input, std::unique_ptr<float[]>& eigens, int dim){
	char JOBZ = 'V';
	char UPLO = 'U';
	int LWORK = 4 * dim;
	unique_ptr<float[]> WORK = make_unique<float[]>(LWORK);
	int INFO;
	ssyev_(&JOBZ, &UPLO, &dim, input.get(), &dim, eigens.get(), WORK.get(), &LWORK, &INFO);
}

template<>
void Eigens_Symm_Wrapper(std::unique_ptr<double[]>& input, std::unique_ptr<double[]>& eigens, int dim){
	char JOBZ = 'V';
	char UPLO = 'U';
	int LWORK = 4 * dim;
	unique_ptr<double[]> WORK = make_unique<double[]>(LWORK);
	int INFO;
	dsyev_(&JOBZ, &UPLO, &dim, input.get(), &dim, eigens.get(), WORK.get(), &LWORK, &INFO);
}
//declarations:
//Eigenvectors
template void Eigens_Wrapper(unique_ptr<double[]>&, unique_ptr<double[]>&, unique_ptr<double[]>&, unique_ptr<double[]>&, int);
template void Eigens_Wrapper(unique_ptr<float[]>&, unique_ptr<float[]>&, unique_ptr<float[]>&, unique_ptr<float[]>&, int);
template void Eigens_Wrapper(unique_ptr<complex<float>[]>&, unique_ptr<complex<float>[]>&, unique_ptr<complex<float>[]>&, unique_ptr<complex<float>[]>&, int);
template void Eigens_Wrapper(unique_ptr<complex<double>[]>&, unique_ptr<complex<double>[]>&, unique_ptr<complex<double>[]>&, unique_ptr<complex<double>[]>&, int);
//Matrix Multiplication
template void Matmul_Wrapper(char, char, unique_ptr<float[]>&, unique_ptr<float[]>&, unique_ptr<float[]>&, int, int, int, float);
template void Matmul_Wrapper(char, char, unique_ptr<double[]>&, unique_ptr<double[]>&, unique_ptr<double[]>&, int, int, int, double);
template void Matmul_Wrapper(char, char, unique_ptr<complex<float>[]>&, unique_ptr<complex<float>[]>&, unique_ptr<complex<float>[]>&, int, int, int, complex<float>);
template void Matmul_Wrapper(char, char, unique_ptr<complex<double>[]>&, unique_ptr<complex<double>[]>&, unique_ptr<complex<double>[]>&, int, int, int, complex<double>);

void Eigens_Symm_Wrapper(std::unique_ptr<float[]>&, std::unique_ptr<float[]>&, int);
void Eigens_Symm_Wrapper(std::unique_ptr<double[]>&, std::unique_ptr<double[]>&, int);
