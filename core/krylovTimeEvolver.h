#pragma once

#include <boost/math/quadrature/tanh_sinh.hpp>

/*
#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t
#include <mkl.h>
#elif defined USE_OPENBLAS
#include <cblas.h>
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#endif
*/
#include <cmath>
#include <complex>
#include <memory>
#include <thread>

#include "matrixDataTypes.h"
#include "krylovObservables.h"
#include "krylovLogger.h"


struct krylovReturn
{
    TE::matrix* sampling;
    std::complex<double>* evolvedState;
    double err;
    size_t n_steps;
	size_t dim;
    size_t krylovDim;
	size_t nSamples;
    int statusCode;

/* Status code has the following meaning: 
1-digit codes mean succes: 0 (everything in order, nothing special happened), 1 (lucky breakdown), 2 (computed analytic error is smaller than estimate of numerical error, which in turn is bigger than requested error; so desired error bound is probably respected) 
more than 1 digit means failure: 10 (computation of error may be  spoiled due to numerical roundoff), 11 (requested tolerance seems unreachable because of roundoff errors), 20 (desired accuracy of numerical integral could not be achieved), 30 (norm of vector deviates significantly from 1), 100 (multiple of these errors)
 */
    krylovReturn(unsigned int nbObservables, unsigned int Hsize, unsigned int nbSamples, int status)
    {
        err = 0; n_steps = 0; krylovDim = 0; dim = Hsize; nSamples = nbSamples; statusCode = status;
        if(nbObservables == 0 && (nSamples * Hsize * sizeof(std::complex<double>) > std::pow(2.,34.)))
        {
            std::cerr << "Requested output would be too large" << std::endl;
            exit(1);
        }
        evolvedState = new std::complex<double>[Hsize];
        if (nbSamples > 0)
        {
            if (nbObservables == 0)
                sampling = new TE::matrix(Hsize, nbSamples);
            else {
                sampling = new TE::matrix(nbObservables, nbSamples);
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
    krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int mm, std::vector<krylovBasicObservable*>  observables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm = true, bool fastIntegration = false, bool progressBar = false);
    krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, std::vector<krylovBasicObservable*>  observables, smatrix* Ham);
    krylovReturn* timeEvolve();
    ~krylovTimeEvolver();

    //sampled values of observables
    TE::matrix* samplings;

    //options
    bool checkNorm, fastIntegration, progressBar;
    std::complex<double> expFactor;
    double tol; size_t m;

    void changeLogLevel(krylovLogger::loggingLevel level);

    
protected:
    int findMaximalStepSize(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double t_step, double t_step_max, int n_s_min, double numericalErrorEstimate, bool increaseStep, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet);
    void sample();
    bool arnoldiAlgorithm(double tolRate, TE::matrix* H, TE::matrix* V, double* h, size_t* m_hbd);
    double integrateError(double a, double b, std::complex<double>* T, std::complex<double>* spectrumH, double h, int method, double tolRate, bool& successful);
    void printProgress(float prog);
    void progressBarThread();


    std::complex<double>* expKrylov(double t, std::complex<double>* T, std::complex<double>* spectrumH);

    
    //Input date
    double t; size_t Hsize;
    double samplingStep; 
    int nbObservables;
    TE::smatrix* Ham;
    std::vector<krylovBasicObservable*>  obsVector;

    //Printing and Logging
    std::thread pBThread;
    krylovLogger logger;
    
    //Determined by input data
    size_t n_samples;
    double matrixNorm;
    double vectorNorm;

    //Internal variables
    boost::math::quadrature::tanh_sinh<double> integ;
    int integrationMethodLong, integrationMethodShort;
    double termination;
    
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
};
