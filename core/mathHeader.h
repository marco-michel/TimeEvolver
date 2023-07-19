#pragma once



#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t
#include <mkl.h>
#include <mkl_spblas.h>
#endif

#ifdef USE_ARMADILLO
#define ARMA_BLAS_CAPITALS
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#include <armadillo>
#endif


#ifdef  USE_OPENBLAS
#include <cblas.h>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#endif
