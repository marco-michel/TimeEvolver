#pragma once


#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

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
    void destroyOptimizeInput();
    void sample();
    void findMaximalStepSize(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double s_0, double t_step_max, int n_s_min, double numericalErrorEstimate, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet);
    bool arnoldiAlgorithm(double tolRate, matrix* H, matrix* V, double* h, size_t* m_hbd);
    
    //Input date
    double t; size_t Hsize;
    double samplingStep; double tol; size_t m;
    smatrix** observables; int nbObservables;
    smatrix* Ham;
    std::complex<double> expFactor;
    bool checkNorm;
    //Determined by input data
    size_t n_samples;
    
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
    size_t index_samples;

    //Useful constants
    static constexpr std::complex<double> one = std::complex<double>(1.0,0.0);
    static constexpr std::complex<double> zero = std::complex<double>(0.0,0.0);
    std::complex<double>* e_1;
};
