#pragma once


#ifdef USE_MKL
#include <mkl.h>
#include <mkl_spblas.h>
#endif

#include "matrixDataTypes.h"

using namespace TE;



enum obsType { VOID_TYPE_OBS, VECTOR_TYPE_OBS, SPARSE_MATRIX_TYPE_OBS, MATRIX_TYPE_OBS };

class krylovBasicObservable
{
public:
    krylovBasicObservable(const std::string& name) : obs_name(name), dim(0), type(VOID_TYPE_OBS) {}
    ~krylovBasicObservable() {}
    virtual std::complex<double> expectation(std::complex<double>* vec, int len) = 0;
    obsType retType();
    std::string retName();

    static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);

protected:
    obsType type;
    std::string obs_name;
    size_t dim;
};



class krylovVectorObservable : public krylovBasicObservable
{
public:
    krylovVectorObservable(const std::string& name, std::complex<double>* obs, size_t len);
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<std::complex<double>[]> obs;
};



class krylovSpMatrixObservable : public krylovBasicObservable
{
public:
    krylovSpMatrixObservable(const std::string& name, smatrix* obs);
    ~krylovSpMatrixObservable();
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<smatrix> obs;
    std::complex<double>* tmpBlasVec;
#ifdef USE_MKL
    sparse_matrix_t* ObsOpt;
    matrix_descr descriptorObs;
#endif
};



class [[deprecated("Dense matrix observables are not fully tested yet. Please use with care.")]] krylovMatrixObservable : public krylovBasicObservable
{
    krylovMatrixObservable(const std::string& name, matrix* obs);
    ~krylovMatrixObservable();
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<matrix> obs;
    std::complex<double>* tmpBlasVec;
};