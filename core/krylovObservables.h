#pragma once

#include <complex>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "mathHeader.h"
#include "matrixDataTypes.h"
#include "parameter.h"


/**
* Exception raised when an observable requests the TimeEvolver to stop at the current timestep.
*/
class requestStopException : public std::exception {
public:
    virtual const char* what() const throw();
    requestStopException() = default;
};


using namespace TE;


/**
* Different types of observable classes describing the underlying data structure
*/
enum obsType { VOID_TYPE_OBS, VECTOR_TYPE_OBS, SPARSE_MATRIX_TYPE_OBS, MATRIX_TYPE_OBS };



/**
* Base observable class used for defining observables, computing expectation values as well as manage output to file
*/
class krylovBasicObservable
{
public:
    krylovBasicObservable(const std::string& name) : obs_name(name), dim(0), numSamples(0), sampleIndex(0), type(VOID_TYPE_OBS), expectationValues(nullptr) {}
    krylovBasicObservable(const std::string& name, std::vector<double> values);
    virtual ~krylovBasicObservable();
    virtual std::complex<double> expectation(std::complex<double>* vec, int len) = 0;
    obsType retType();
    std::string retName();
    double* retexpectationValues();
    size_t retnumSamples();
    virtual void initializeResultArray(size_t size);
    size_t resetResultArray();

    static void saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>> &obs_list, parameter_list& para, const std::string& name);


    static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);

protected:
    std::string obs_name;
    size_t dim;
    size_t numSamples;
    size_t sampleIndex;
    obsType type;
    double* expectationValues;
};

/**
* Derived observable class only used for easier file output. Can't be used for sampling in the TimeEvolver. Calling the expectation function will result in a runtime error. 
*/
class krylovOutputObservable : public krylovBasicObservable {
public:
    krylovOutputObservable(const std::string& name, std::vector<double> values);
    std::complex<double> expectation(std::complex<double>* vec, int len);
};

/**
* Derviced observable class. The observable is represented as a vector and the expecation value is computed as a dot product: <Obs|state>
*/
class krylovVectorObservable : public krylovBasicObservable
{
public:
    krylovVectorObservable(const std::string& name, std::complex<double>* obs, size_t len);
    ~krylovVectorObservable() {}
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<std::complex<double>[]> obs;
};


/**
* Derived observable class. The observable is represented as a sparse matrix and the expectation value is computed as <Obs|state|Obs>
*/
class krylovSpMatrixObservable : public krylovBasicObservable
{
public:
    krylovSpMatrixObservable(const std::string& name, std::unique_ptr<smatrix> obs);
    ~krylovSpMatrixObservable();
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<smatrix> obs;
    std::complex<double>* tmpBlasVec;
};


/**
* Derived observable class. The observable is represented as a (dense) matrix and the expectation value is computed as <Obs|state|Obs>
*/
class krylovMatrixObservable : public krylovBasicObservable
{
    krylovMatrixObservable(const std::string& name, std::unique_ptr<matrix> obs);
    ~krylovMatrixObservable();
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<matrix> obs;
    std::complex<double>* tmpBlasVec;
};