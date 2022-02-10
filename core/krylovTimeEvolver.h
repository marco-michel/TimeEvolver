#pragma once


#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <boost/math/quadrature/tanh_sinh.hpp>

#include <mkl.h>
#include <mkl_spblas.h>

#include <cmath>
#include <complex>
#include "matrixDataTypes.h"

struct krylovReturn
{
    matrix* sampling;
    std::complex<double>* evolvedState;
    double err;
    size_t n_steps;
	size_t dim;
	size_t nSamples;

    krylovReturn(unsigned int nbObservables, unsigned int Hsize, unsigned int nbSamples)
    {
        err = 0; n_steps = 0; dim = Hsize; nSamples = nbSamples;
        if(nbObservables == 0 && (nSamples * Hsize * sizeof(std::complex<double>) > std::pow(2.,34.)))
        {
            std::cerr << "Requested output would be too large" << std::endl;
            exit(1);
        }
        evolvedState = new std::complex<double>[Hsize];
        if(nbSamples > 0)
        {
            if(nbObservables == 0)
                sampling = new matrix(Hsize, nbSamples);
            else {
                sampling = new matrix(nbObservables, nbSamples);
            }
        }
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
    krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int mm, smatrix** observables, int nbObservables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm);
    krylovReturn* timeEvolve();
    ~krylovTimeEvolver();

    //sampled values of observables
    matrix* samplings;
    
protected:
    void optimizeInput();
    void findMaximalStepSize2(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double t_step, double t_step_max, int n_s_min, double numericalErrorEstimate, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet, bool increaseStep);
    void destroyOptimizeInput();
    void sample();
    void findMaximalStepSize(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double s_0, double t_step_max, int n_s_min, double numericalErrorEstimate, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet);
    bool arnoldiAlgorithm(double tolRate, matrix* H, matrix* V, double* h, size_t* m_hbd);
    double integrateError(double a, double b, std::complex<double>* T, std::complex<double>* spectrumH, double h);

    double errorKernel(double t, std::complex<double>* T, std::complex<double>* spectrumH, double h);

    std::complex<double>* expKrylov(double t, std::complex<double>* T, std::complex<double>* spectrumH);

    
    //Input date
    double t; size_t Hsize;
    double samplingStep; double tol; size_t m;
    smatrix** observables; int nbObservables;
    smatrix* Ham;
    std::complex<double> expFactor;
    bool checkNorm;
    
    //Determined by input data
    size_t n_samples;
    double matrixNorm;

    //Internal variables
    boost::math::quadrature::tanh_sinh<double> integ;
    
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
    std::complex<double>* tmpintKernelT;
    size_t index_samples;

    //Useful constants
    static constexpr std::complex<double> one = std::complex<double>(1.0,0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0,0.0);
    std::complex<double>* e_1;
};
