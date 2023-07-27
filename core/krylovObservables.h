#pragma once

#include <complex>
#include <fstream>
#include <iostream>
#include <memory>

#include "mathHeader.h"
#include "matrixDataTypes.h"
#include "parameter.h"






using namespace TE;


enum obsType { VOID_TYPE_OBS, VECTOR_TYPE_OBS, SPARSE_MATRIX_TYPE_OBS, MATRIX_TYPE_OBS };

class krylovBasicObservable
{
public:
    krylovBasicObservable(const std::string& name) : obs_name(name), dim(0), type(VOID_TYPE_OBS), numSamples(0), sampleIndex(0), expectationValues(nullptr) {}
    virtual ~krylovBasicObservable();
    virtual std::complex<double> expectation(std::complex<double>* vec, int len) = 0;
    obsType retType();
    std::string retName();
    double* retexpectationValues();
    void initializeResultArray(size_t size);

    static void saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>> &obs_list, parameter_list& para, const std::string& name);


    static constexpr std::complex<double> one = std::complex<double>(1.0, 0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0, 0.0);


    std::string obs_name;
    size_t dim;
    size_t numSamples;



protected:
    size_t sampleIndex;
    obsType type;
    double* expectationValues;
};



class krylovVectorObservable : public krylovBasicObservable
{
public:
    krylovVectorObservable(const std::string& name, std::complex<double>* obs, size_t len);
    ~krylovVectorObservable() {}
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<std::complex<double>[]> obs;
};



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



class krylovMatrixObservable : public krylovBasicObservable
{
    krylovMatrixObservable(const std::string& name, std::unique_ptr<matrix> obs);
    ~krylovMatrixObservable();
    std::complex<double> expectation(std::complex<double>* vec, int len);

private:
    std::unique_ptr<matrix> obs;
    std::complex<double>* tmpBlasVec;
};