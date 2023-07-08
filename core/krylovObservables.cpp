
#include <iostream>
#include <complex>
#include <memory>

#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t
#include <mkl.h>
#elif defined USE_OPENBLAS
#include <cblas.h>
//#include <lapacke.h>
#endif

#include "matrixDataTypes.h"
#include "krylovObservables.h"

using namespace TE;


std::string krylovBasicObservable::retName()
{
	return obs_name;
}

obsType krylovBasicObservable::retType()
{
	return type;
}


krylovMatrixObservable::krylovMatrixObservable(const std::string& name, matrix* obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = MATRIX_TYPE_OBS;
	if (dim > 0)
		tmpBlasVec = new std::complex<double>[dim];

	obs = std::make_unique<matrix>(*obser);

}

krylovMatrixObservable::~krylovMatrixObservable()
{
	if (dim > 0)
		delete[] tmpBlasVec;
}

std::complex<double> krylovMatrixObservable::expectation(std::complex<double>* vec, int len) //requires testing
{
	if (len != dim)
	{
		std::cerr << "Incompatible dimensions" << std::endl;
		exit(1);
	}
	std::complex<double> observall;
	cblas_zgemv(CblasColMajor, CblasNoTrans, dim, dim, &one, obs->values, dim, vec, 1, &zero, tmpBlasVec, 1);
	cblas_zdotc_sub(len, vec, 1, tmpBlasVec, 1, &observall);
	return observall;
}

krylovSpMatrixObservable::krylovSpMatrixObservable(const std::string& name, smatrix* obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = SPARSE_MATRIX_TYPE_OBS;
	obs = std::make_unique<smatrix>(*obser); //make owned copy (I don't know if that's the best way)
	obs->initialize();

	tmpBlasVec = new std::complex<double>[dim];
}

krylovSpMatrixObservable::~krylovSpMatrixObservable()
{
	if (dim > 0)
	{
		delete[] tmpBlasVec;
	}
}

std::complex<double> krylovSpMatrixObservable::expectation(std::complex<double>* vec, int len)
{

	if (len != dim)
	{
		std::cerr << "Incompatible dimensions" << std::endl;
		exit(1);
	}

	std::complex<double> observall;
	obs->spMV(one, vec, tmpBlasVec);
	cblas_zdotc_sub(len, vec, 1, tmpBlasVec, 1, &observall);

	return observall;

}

//this would not be necessary with C++17

constexpr std::complex<double> krylovBasicObservable::zero;
constexpr std::complex<double> krylovBasicObservable::one;

krylovVectorObservable::krylovVectorObservable(const std::string& name, std::complex<double>* obser, size_t len) : krylovBasicObservable(name)
{
	dim = len;
	type = VECTOR_TYPE_OBS;
	obs = std::make_unique<std::complex<double>[]>(len);
	cblas_zcopy(dim, obser, 1, obs.get(), 1);
}

std::complex<double> krylovVectorObservable::expectation(std::complex<double>* vec, int len)
{
	if (len != dim)
	{
		std::cerr << "Incompatible dimensions" << std::endl;
		exit(1);
	}
	std::complex<double> observall;
	cblas_zdotc_sub(len, vec, 1, obs.get(), 1, &observall);
	std::complex<double> observallreturn;
	observallreturn.imag(0);
	observallreturn.real(std::norm(observall));
	return observallreturn; //returns the squared magnitued
}
