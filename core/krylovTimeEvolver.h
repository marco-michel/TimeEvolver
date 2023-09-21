#pragma once

#include <complex>
#include <limits>
#include <algorithm>
#include <thread>
#include <sstream>
#include <vector>
#include <memory>
#include <exception>
#include <stdexcept>
#include <iomanip>

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp> 
#include <boost/throw_exception.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include "mathHeader.h"
#include "matrixDataTypes.h"
#include "krylovObservables.h"
#include "krylovLogger.h"




/**
* Return type of krylovTimeEvolver containing general information as well as the original Hamiltonian and the observables containing the actual expectation values.
* Note that the matrix representations might have been changed in the course of the time evolution 
*/
class krylovReturn
{
public:
    std::complex<double>* evolvedState;
    double err; double evolvedTime;
    size_t numSamples;
    size_t n_steps;
	size_t dim;
    size_t krylovDim;
    /* Status code has the following meaning:
    1-digit codes mean succes: 0 (everything in order, nothing special happened), 1 (lucky breakdown), 2 (computed analytic error is smaller than estimate of numerical error, which in turn is bigger than requested error; so desired error bound is probably respected)
    more than 1 digit means failure: 10 (computation of error may be  spoiled due to numerical roundoff), 11 (requested tolerance seems unreachable because of roundoff errors), 20 (desired accuracy of numerical integral could not be achieved), 30 (norm of vector deviates significantly from 1), 100 (multiple of these errors)
    */
    int statusCode;
    std::vector<std::unique_ptr<krylovBasicObservable>> observableList;
    std::unique_ptr<smatrix> hamiltonianMatrix;

    krylovReturn(unsigned int Hsize, int status);
    ~krylovReturn();
};

/**
* Main class responsible for computing the actual time evolution.
*/
class krylovTimeEvolver
{
public:
    krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, double tol, int mm, std::vector<std::unique_ptr<krylovBasicObservable>>  observables, std::unique_ptr<smatrix> HamIn, double expFactor, bool fastIntegration, bool progressBar);
    krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, std::vector<std::unique_ptr<krylovBasicObservable>> observables, std::unique_ptr<smatrix> Ham);
    krylovReturn* timeEvolve();
    ~krylovTimeEvolver();

    //sampled values of observables


    //options
    bool fastIntegration, progressBar;
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
    krylovReturn* generateReturn();


    std::complex<double>* expKrylov(double t, std::complex<double>* T, std::complex<double>* spectrumH);

    
    //Input date
    double t; size_t Hsize;
    double samplingStep; 
    int nbObservables;
    std::unique_ptr<smatrix> Ham;
    std::vector<std::unique_ptr<krylovBasicObservable>>  obsVector;

    //Printing and Logging
    std::thread pBThread;
    std::atomic<bool> stop_printing;
    krylovLogger logger;
    
    //Determined by input data
    size_t n_samples;
    double matrixNorm;

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
    //Time of current state
    double t_now = 0.;
    //Number of steps so far
    int n_steps = 0;
    //Latest time at which sampling has happened so far
    double t_sampling = 0.;
    //Total accumulated error so far
    double err = 0.;
    //Track if an error occured within timeEvolve(). A value unequal 0 indicates that numerical result likely violates error bound.
    int statusCode = 0;

    //Useful constants
    static constexpr std::complex<double> one = std::complex<double>(1.0,0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0,0.0);
    std::complex<double>* e_1;
    static const int pBarWidth = 70;
};
