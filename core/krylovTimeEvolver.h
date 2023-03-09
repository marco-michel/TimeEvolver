#pragma once

#include <boost/math/quadrature/tanh_sinh.hpp>

#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <mkl.h>
#include <mkl_spblas.h>

#include <cmath>
#include <complex>
#include <memory>

#include "matrixDataTypes.h"


enum obsType {VOID_TYPE_OBS, VECTOR_TYPE_OBS, SPARSE_MATRIX_TYPE_OBS, MATRIX_TYPE_OBS };

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
    matrix_descr descriptorObs;
    std::complex<double>* tmpBlasVec;
    sparse_matrix_t* ObsOpt;
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

struct krylovReturn
{
    matrix* sampling;
    std::complex<double>* evolvedState;
    double err;
    size_t n_steps;
	size_t dim;
	size_t nSamples;
    int statusCode;

/* Status code has the following meaning: 
1-digit codes mean succes: 0 (everything in order, nothing special happened), 1 (lucky breakdown), 2 (computed analytic error is smaller than estimate of numerical error, which in turn is bigger than requested error; so desired error bound is probably respected) 
more than 1 digit means failure: 10 (computation of error may be  spoiled due to numerical roundoff), 11 (requested tolerance seems unreachable because of roundoff errors), 20 (desired accuracy of numerical integral could not be achieved), 30 (norm of vector deviates significantly from 1), 100 (multiple of these errors)
 */
    krylovReturn(unsigned int nbObservables, unsigned int Hsize, unsigned int nbSamples, int status)
    {
        err = 0; n_steps = 0; dim = Hsize; nSamples = nbSamples; statusCode = status;
        if(nbObservables == 0 && (nSamples * Hsize * sizeof(std::complex<double>) > std::pow(2.,34.)))
        {
            std::cerr << "Requested output would be too large" << std::endl;
            exit(1);
        }
        evolvedState = new std::complex<double>[Hsize];
        if (nbSamples > 0)
        {
            if (nbObservables == 0)
                sampling = new matrix(Hsize, nbSamples);
            else {
                sampling = new matrix(nbObservables, nbSamples);
            }
        }
        else
            sampling = nullptr;
    }

	~krylovReturn()
	{
            delete[] evolvedState;
            if (nSamples > 0)
				delete sampling;
	}
};


class krylovTimeEvolver
{
public:
    krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int mm, smatrix** observables, int nbObservables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm= true, bool fastIntegration = false);
    krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int mm, std::vector<std::unique_ptr<krylovBasicObservable>>  observables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm = true, bool fastIntegration = false, bool progressBar = false);
    krylovReturn* timeEvolve();
    ~krylovTimeEvolver();

    //sampled values of observables
    matrix* samplings;
    
protected:
    void optimizeInput();
    int findMaximalStepSize(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double t_step, double t_step_max, int n_s_min, double numericalErrorEstimate, bool increaseStep, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet);
    void destroyOptimizeInput();
    void sample();
    bool arnoldiAlgorithm(double tolRate, matrix* H, matrix* V, double* h, size_t* m_hbd);
    double integrateError(double a, double b, std::complex<double>* T, std::complex<double>* spectrumH, double h, int method, double tolRate, bool& successful);
    void printProgress(float prog);

    std::complex<double>* expKrylov(double t, std::complex<double>* T, std::complex<double>* spectrumH);

    
    //Input date
    double t; size_t Hsize;
    double samplingStep; double tol; size_t m;
    smatrix** observables; int nbObservables;
    smatrix* Ham;
    std::complex<double> expFactor;
    bool checkNorm, fastIntegration, progressBar;
    
    //Determined by input data
    size_t n_samples;
    double matrixNorm;

    //Internal variables
    boost::math::quadrature::tanh_sinh<double> integ;
    int integrationMethodLong, integrationMethodShort;
    double termination;
    
    //variables for mkl-library
    sparse_matrix_t** ObsOpt;
    sparse_matrix_t* HamOpt;
    matrix_descr descriptor;
    matrix_descr descriptorObs;
    
    //temporary variables shared by different functions
    std::complex<double>* currentVec;
    std::complex<double>* sampledState;
    std::complex<double>* tmpBlasVec;
    std::complex<double>* tmpKrylovVec1;
    std::complex<double>* tmpKrylovVec2;
    std::complex<double>* tmpintKernelExp;
    std::complex<double>* tmpintKernelExp1;
    std::complex<double>* tmpintKernelExp2;
    std::complex<double>* tmpintKernelExp3;
    std::complex<double>* tmpintKernelT;
    size_t index_samples;

    //Useful constants
    static constexpr std::complex<double> one = std::complex<double>(1.0,0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0,0.0);
    std::complex<double>* e_1;
    static const int pBarWidth = 70;

    //New observables
    bool obsComputeExpectation;
    std::vector<std::unique_ptr<krylovBasicObservable>>  obsVector;
};
