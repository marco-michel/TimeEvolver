/* Brief summary :
 For given Hamiltonian, time evolve a given initial state, while also saving results at intermediate time steps.
 This is done in such a way that at each step, the norm difference of the state obtained by numerical time evolution deviates from the true result at most by a prescribed error bound.
 The numerical routine is based on iteratively finding a Krylov-subspaces, on which time evolution can be perform for a small time step.
 The time-evolved vector then becomes the initial vector for constructing the subsequent Krylov-subspace.
 */

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss.hpp> 
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <complex.h>
#include <time.h>
#include <complex>
#include <limits>
#include <algorithm>
#include <iostream>

#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <mkl.h>
#include <mkl_spblas.h>

#include "matrixDataTypes.h"
#include "krylovTimeEvolver.h"

/**
 * All data required for the numerical time evolution is needed.
 * @param t The time interval over which the state should be time evolved.
 * @param Hsize The size of the full Hilbert space
 * @param v The initial state that should be time evolved
 * @param samplingStep The time interval after which the values of the observables should be determined
 * @param tol The maximal admissible error (norm difference between result of numerical and true time evolution)
 * @param m The size of the Krylov subspaces
 * @param observables The observables that are to be sampled
 * @param nbObservables The number of observables
 * @param Ham The full Hamiltonian
 * @param expFactor The scalar factor multiplying the Hamiltonian in the time evolution (usually -i)
 * @param checkNorm Whether or not it should be check that the time evolved state has unit norm
 */
krylovTimeEvolver::krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int m, smatrix** observables, int nbObservables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm)
{
    if(Hsize == 0)
    {
        std::cerr << "Invalid Hilbertspace dimension" << std::endl;
        exit(1);
    }
    
    this->t = t; this->Hsize = Hsize; this->samplingStep = samplingStep; this->tol = tol; this->m  = std::min<size_t>(m, Hsize);
    this->observables = observables; this->nbObservables = nbObservables; this->Ham = Ham; this->expFactor = expFactor; this->checkNorm = checkNorm;
    ObsOpt = nullptr; HamOpt = nullptr;

	matrixNorm = Ham->normInf();
    
    //number of sampling steps
    n_samples = (size_t) floor(t / samplingStep) + 1;
    
    if (nbObservables == 0)
        samplings = new matrix(Hsize, n_samples);
    else
        samplings = new matrix(nbObservables, n_samples);
    
    //The state at the current time
    currentVec = new std::complex<double>[Hsize];
    cblas_zcopy(Hsize, v, 1, currentVec, 1);
    //The state to be sampled
    sampledState = new std::complex<double>[Hsize];
    cblas_zcopy(Hsize, v, 1, sampledState, 1);
    //A temporary vector of size Hsize
    tmpBlasVec = new std::complex<double>[Hsize];
    //Two temporary vectors of size m
    tmpKrylovVec1 = new std::complex<double>[m];
    tmpKrylovVec2 = new std::complex<double>[m];
	tmpintKernelExp = new std::complex<double>[m];
	tmpintKernelExp1 = new std::complex<double>[m];
	tmpintKernelExp2 = new std::complex<double>[m];
	tmpintKernelT = new std::complex<double>[m];

    descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
    descriptorObs.diag = SPARSE_DIAG_NON_UNIT;
    descriptorObs.type = SPARSE_MATRIX_TYPE_GENERAL;

    index_samples = 0;
    e_1 = new std::complex<double>[m];
    e_1[0].real(1);

}

krylovTimeEvolver::~krylovTimeEvolver()
{
    delete[] tmpBlasVec;
    delete[] tmpKrylovVec1;
    delete[] tmpKrylovVec2;
    delete[] sampledState;
    delete[] currentVec;
	delete[] tmpintKernelExp;
	delete[] tmpintKernelExp1;
	delete[] tmpintKernelExp2;
	delete[] tmpintKernelT;
    delete samplings;
    if(HamOpt != nullptr)
        delete HamOpt;
    if(HamOpt != nullptr)
    {
        for (int i = 0; i != nbObservables; i++)
            delete ObsOpt[i];
        delete[] ObsOpt;
    }
    delete[] e_1;
}

constexpr std::complex<double> krylovTimeEvolver::one;
constexpr std::complex<double> krylovTimeEvolver::zero;

/**
 * Computes and saves values of observables for current state ('sampledState')
 */
void krylovTimeEvolver::sample()
{
    if (nbObservables == 0)
        cblas_zcopy(Hsize, sampledState, 1, samplings->values + index_samples*Hsize, 1);
    else
    {
        std::complex<double> observall;
        for (int i = 0; i < nbObservables; i++)
        {
            mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,one,*ObsOpt[i],descriptorObs,sampledState,zero,tmpBlasVec);
            cblas_zdotc_sub(Hsize, sampledState, 1, tmpBlasVec, 1, &observall);
            *(samplings->values + i + index_samples * nbObservables) = observall.real();
        }
    }
    index_samples++;
}

/**
 * Destroys variables created in krylovTimeEvolver::optimizeInput
 */
void krylovTimeEvolver::destroyOptimizeInput()
{
    if(HamOpt != nullptr)
        mkl_sparse_destroy(*HamOpt);
    if(nbObservables > 0)
    {
        for(int i = 0; i != nbObservables; i++)
            mkl_sparse_destroy(*ObsOpt[i]);
    }
}

/**
 * Brings input parameters in form needed for mkl-routines that perform operations on large matrices and vectors
 */
void krylovTimeEvolver::optimizeInput()
{
    if (Ham == nullptr)
        return;
    
    sparse_status_t mklStatus;
    matrix_descr type; type.type = SPARSE_MATRIX_TYPE_GENERAL;
    
    HamOpt = new sparse_matrix_t;
    ObsOpt = new sparse_matrix_t*[nbObservables];
    for (int i = 0; i != nbObservables; i++)
        ObsOpt[i] = new sparse_matrix_t;
    
    if (Ham->numValues != 0)
        mklStatus = mkl_sparse_z_create_coo(HamOpt, SPARSE_INDEX_BASE_ZERO, Ham->m, Ham->n, Ham->numValues, Ham->rowIndex, Ham->columns, Ham->values);
    
    mklStatus = mkl_sparse_convert_csr(*HamOpt, SPARSE_OPERATION_NON_TRANSPOSE, HamOpt);
    mklStatus = mkl_sparse_order(*HamOpt);
    mklStatus = mkl_sparse_set_mv_hint(*HamOpt, SPARSE_OPERATION_NON_TRANSPOSE, type, (size_t)std::llabs(std::llround(m*t)));
    mklStatus = mkl_sparse_set_memory_hint(*HamOpt, SPARSE_MEMORY_AGGRESSIVE);
    mklStatus = mkl_sparse_optimize(*HamOpt);
    
    for (int i = 0; i != nbObservables; i++)
    {
        mklStatus = mkl_sparse_z_create_coo(ObsOpt[i], SPARSE_INDEX_BASE_ZERO, observables[i]->m, observables[i]->n, observables[i]->numValues, observables[i]->rowIndex, observables[i]->columns, observables[i]->values);
        if(mklStatus != SPARSE_STATUS_SUCCESS)
        {
            std::cerr << "Could not process Matrix representation of observable " << i << ". Empty matrices can not be processed" << std::endl;
            exit(1);
        }
    }
}


void krylovTimeEvolver::findMaximalStepSize2(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double t_stepSuggestion, double t_step_max, int n_s_min, double numericalErrorEstimate, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet, bool increaseStep)
{
	//Maximal number of substepreductions to meet tolerance
	unsigned int GO_MAX = 100;
	unsigned int nbReductions = 0;

	double t_step = t_stepSuggestion;

	//Increasing step size. Used in the first go through when there is no good guess for the step size. 
	if (increaseStep)
	{
		while (integrateError(0, 2*t_step, T, spectrumH, h) < tolRate * t_step && t_step < t_step_max)
			t_step *= 2.0;
	}

	double err_step = integrateError(0, t_step, T, spectrumH, h);
	
	while (err_step > tolRate * t_step)
	{
		nbReductions++;
		t_step = t_step / 2.0;
		err_step = integrateError(0, t_step, T, spectrumH, h);
		if (nbReductions == GO_MAX)
		{
			std::cerr << "Error: No small enough time step found to meet tolerance requirements." << std::endl;
			exit(-1);
		}
	}

	double s = t_step / n_s_min;
	int n_s = 0;
	double deltaError = 0;

	//Maximal number of steps to reach t_step_max
	double n_s_max = ceil(t_step_max / s);
	if (n_s_max < n_s_min) {
		n_s_max = n_s_min + 1;
		s = t_step_max / n_s_max;
	}

	while (err_step + deltaError < tolRate * (t_step + n_s * s) && n_s <= n_s_max - 1)
	{
		deltaError += integrateError(t_step + n_s * s, t_step + (n_s + 1) * s, T, spectrumH, h);
		n_s++;
	}

	t_step += (n_s - 1) * s;
	err_step += deltaError;

	if (err_step < numericalErrorEstimate && (n_s != n_s_max)) {
		std::cerr
			<< "CRITICAL WARNING: the computed error bound "
			<< err_step
			<< "was smaller than the estimate of the numerical error "
			<< numericalErrorEstimate
			<< ". THE DESIRED ERROR BOUND WILL LIKELY BE VIOLATED. (Remaining time "
			<< t_step_max
			<< ") Restart with bigger error bound or smaller time."
			<< std::endl;
		err_step = numericalErrorEstimate;
	}

	*t_stepRet = t_step;
	*err_stepRet = err_step;

	std::complex<double>* w_Krylov = expKrylov(t_step, T, spectrumH);
	cblas_zcopy(m, w_Krylov, 1, w_KrylovRet, 1);

}


/**
 * Given a Krylov subspace, this function finds out how far (i.e. for what time step) one can use it without exceeding the prescribed error bound. To this end, it uses small substeps and computes an error bound for each substep. This function only operates in the Krylov subspace, i.e. with vectors and matrices of dimension m
 * @param T The orthogonal transformation matrix for the Hamiltonian in the Krylov subspace
 * @param spectrumH The eigenvalues of the Hamiltonian in the Krylov subspace
 * @param h The last element of the Arnoldi algorithm (it's needed for the computation of the error)
 * @param tolRate The maximal admissible error rate (i.e. error per time)
 * @param s The initial substep size
 * @param t_step_max The maximal stepsize (in case tolRate is not exceeded)
 * @param n_s_min The minimal number of substeps (the numerical integration of the error bound gets bad if there are too few substeps)
 * @param numericalErrorEstimate The numerical error associated to the Krylov subspace (this sets a lower bound for the error)
 * @param t_stepRet Returns the time step
 * @param w_KrylovRet Returns the Krylov vector after the time step
 * @param err_stepRet Returns the error of the time step
 */
void krylovTimeEvolver::findMaximalStepSize(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double s, double t_step_max, int n_s_min, double numericalErrorEstimate, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet)
{

	//CONSTANTS
	//Maximal number of substepreductions to meet tolerance
	unsigned int GO_MAX = 100;
	// If the local error rate is bigger than tolRate by this fraction, it is
	// regarded as non-zero and sign changes in the error rate will be monitored
	double SIGNIFICANT_ERROR_FRACTION = 1. / 100;
	// If the number of monitored sign changes of the derivative of the error is
	//bigger than the number of steps by this fraction, a warning will be issued
	double SIGNIFICANT_SIGN_CHANGE_FRACTION = 1. / 7;

	//TEMPORARY VARIABLES
	double t_step = 0;
	bool substepReduction;
	unsigned int nbSignChanges;
	int n_s;
	std::complex<double> *w_Krylov = new std::complex<double>[m];
	std::complex<double> *w_Krylov_Previous = new std::complex<double>[m];
	std::complex<double> *expSpectrum = new std::complex<double>[m];
	std::complex<double> *spectrumHTime = new std::complex<double>[m];

	//Maximal number of steps to reach t_step_max
	double n_s_max = ceil(t_step_max / s);
	if (n_s_max < n_s_min) {
		n_s_max = n_s_min + 1;
		s = t_step_max / n_s_max;
	}

	//In each loop, one tries to reach the minimal number of required substeps n_s_min with a given substepsize s
	//If it is not successful, s is halved
	for (unsigned int go = 1; go < GO_MAX; go++) {
		n_s = 0;
		substepReduction = false;
		double err_step = 0;
		double err_step_Previous = 0;
		double errRate_Previous = 0;
		double errRate = 0;
		nbSignChanges = 0;
		int theSignPrevious = 1;

		//Consecutive substeps are performed until tolRate is exceeded
		for (; n_s <= n_s_max - 1; n_s++) {
			//Compute time-evolved Krylov vector: w_Krylov = T e^((n_s + 1)*s*spectrumH)*T'*e_1
			cblas_zcopy(m, spectrumH, 1, spectrumHTime, 1);
			cblas_zdscal(m, (n_s + 1) * s, spectrumHTime, 1);
			vzExp(m, spectrumHTime, expSpectrum);
			cblas_zgemv(CblasColMajor, CblasConjTrans, m, m, &one, T, m, e_1, 1,
					&zero, tmpKrylovVec1, 1);
			for (unsigned int i = 0; i != m; i++)
				tmpKrylovVec1[i] = tmpKrylovVec1[i] * expSpectrum[i];
			cblas_zgemv(CblasColMajor, CblasNoTrans, m, m, &one, T, m,
					tmpKrylovVec1, 1, &zero, w_Krylov, 1);

			//Compute local error via numerical integration: err_step += s*h*max(|w_Krylov[m - 1]|,|w_Krylov:Previous[m - 1]|)
			errRate = h * std::abs(w_Krylov[m - 1]);
			err_step += s * std::max(errRate, errRate_Previous);
			//To check quality of numerical integration: count how many times errorRate changes from decreasing to increasing
			int theSign = (int) std::signbit(errRate - errRate_Previous);
			if (errRate > SIGNIFICANT_ERROR_FRACTION * tolRate
					&& theSign != theSignPrevious)
				nbSignChanges++;

			t_step = (n_s + 1) * s;

			//Stop loop over substeps if tolRate is exceeded
			if (err_step > (n_s + 1) * s * tolRate) {
				//If there were not enough substeps, try again with smaller substepsize
				if (n_s < n_s_min) {
					substepReduction = true;
					s = fmin(s / 2.0, s * ((n_s + 1.0) / n_s_min));
					n_s_max = ceil(t_step_max / s);
				}
				//If there were enough substeps, the function can terminate
				//The current substep is discarded (since it exceeded tolRate) and the previous substep is returned
				else {
					t_step = n_s * s;
					err_step = err_step_Previous;
					cblas_zcopy(m, w_Krylov_Previous, 1, w_Krylov, 1);
					if (err_step < numericalErrorEstimate) {
						std::cerr
								<< "CRITICAL WARNING: the computed error bound "
								<< err_step
								<< "was smaller than the estimate of the numerical error "
								<< numericalErrorEstimate
								<< ". THE DESIRED ERROR BOUND WILL LIKELY BE VIOLATED. (Remaining time "
								<< t_step_max
								<< ") Restart with bigger error bound or smaller time."
								<< std::endl;
						err_step = numericalErrorEstimate;
					}
				}
				break;
			}

			cblas_zcopy(m, w_Krylov, 1, w_Krylov_Previous, 1);
			errRate_Previous = errRate;
			err_step_Previous = err_step;
			theSignPrevious = theSign;

		}//end n_s-loop

	    //If enough substeps were achieved, the go-loop can be terminated
		if (substepReduction == false) {
			*t_stepRet = t_step;
			*err_stepRet = err_step;
			cblas_zcopy(m, w_Krylov, 1, w_KrylovRet, 1);
			if (nbSignChanges > SIGNIFICANT_SIGN_CHANGE_FRACTION * n_s) {
				std::cout
						<< "WARNING: The sampling may not have been sufficiently dense. "
						<< "The error rate has changed between increasing and decreasing "
						<< nbSignChanges << " times with only " << n_s
						<< " total substeps. (Remaining time " << t_step_max
						<< ")" << std::endl;
			}
			break;
		}

	}//end go-loop

	//If the go-loop was terminated without ever achieving enough substeps, an error is returned
	if (substepReduction == true) {
		std::cerr << "ERROR : Substep " << s
				<< " still too big at remaining time " << t_step_max
				<< " after maximum number " << GO_MAX
				<< " of substep-reuctions reached, without the error staying below the"
				<< " required tolerance" << std::endl;
		exit(1);
	}

	delete[] w_Krylov_Previous;
	delete[] w_Krylov;
	delete[] spectrumHTime;
	delete[] expSpectrum;
}

/**
 * Main rountine: performs time evolution according to the data provided in the constructor
 * @return Returns in particular the samplings for each observable, the state after time evolution and an upper bound on the error
 */
krylovReturn* krylovTimeEvolver::timeEvolve()
{
	//Constants
	//Minimal number of substeps per time step (needed for accuracy of computation of error)
	int N_SUBSTEPS_MIN = 50;
	//After each time step, the optimal step size is computed. To avoid substep reduction because of a too small number of substeps, the optimal step size is multiplied by this number
	double INITIAL_STEP_FRACTION = 0.97;

	optimizeInput();
	if (checkNorm) {
		if (cblas_dznrm2(Hsize, currentVec, 1) - 1.0 > tol) {
			std::cerr << "Norm error in initial vector" << std::endl;
			exit(1);
		}
	}
    
    //Time of current state
	double t_now = 0.;
	//Number of steps so far
	int n_steps = 0;
	//Estimate for initial step size
	double t_step = 1.0 / 10.0;
	//Latest time at which sampling has happened so far
	double t_sampling = 0.;

	//Admissible error rate
	double tolRate = tol / t;
	//Total accumulated error so far
	double err = 0.;
	//Estimate of error due to limited precision of numerical operations
	double numericalErrorEstimate = Hsize * matrixNorm * std::numeric_limits<double>::epsilon();

    //Flag indicating if a happy breakdown has occured
    bool dummy_hbd = false;
    //In case of lucky breakdown, size of Krylov space
    size_t m_hbd;

    //Hessenberg matrix
	matrix *H = new matrix(m, m);
	//Corresponding transformation matrix
	matrix *V = new matrix(Hsize, m);
	//The (m+1,m) element of Hessenberg matrix (needed for computation of error)
	double h = 0;
	//Eigenvalues of Hessenberg matrix
	std::complex<double> *eigenvalues = new std::complex<double>[m];
	//Eigenvectors of Hessenberg matrix
	std::complex<double> *schurvector = new std::complex<double>[m * m];

	//Record observables for initial state
	sample();

    //Main loop
    while (index_samples < n_samples)
    {
		n_steps++;
		t_step = fmin(t - t_now, t_step);
		double err_step = 0;

		//STEP 1: Construct Krylov subspace using Arnoldi algorithm
		dummy_hbd = arnoldiAlgorithm(tolRate, H, V, &h, &m_hbd);

		//Some special adjustments in case of a lucky breakdown, i.e. when projection in Krylov-subspace of dimension m_hbd <= m is exact (within numerical uncertainty)
		//In particular, the time step of the current Krylov space can be arbitarily large in this case
		if (dummy_hbd) {
			t_step = t - t_now;
			matrix *Htmp = new matrix(m_hbd, m_hbd);
			matrix *Vtmp = new matrix(Hsize, m_hbd);
			for (size_t ll = 0; ll != m_hbd * m_hbd; ll++) {
				Htmp->values[ll] = H->values[ll
						+ (m - m_hbd) * (int) std::floor(ll / m_hbd)];
			}
			for (size_t ll = 0; ll != m_hbd * Hsize; ll++) {
				Vtmp->values[ll] = V->values[ll];
			}
			delete H;
			delete V;
			H = Htmp;
			V = Vtmp;
			m = m_hbd;
			err_step = std::max(h,numericalErrorEstimate)*t_step;
		}
		//Finally diagonalize Hessenberg matrix H (since it will be exponentiated many times)
		int infocheck = LAPACKE_zhseqr(LAPACK_COL_MAJOR, 'S', 'I', m, 1, m,
				H->values, m, eigenvalues, schurvector, m);
		if (infocheck != 0) {
			std::cout << "LAPACK Error " << infocheck << std::endl;
		}
		//END STEP 1
        
        //STEP 2: find maximal time step ('t_step') for which current Krylov subspace can be used without violating error bound
		if (dummy_hbd == false)
        {
            //double s_0 = INITIAL_STEP_FRACTION*t_step / N_SUBSTEPS_MIN;
            //findMaximalStepSize(schurvector, eigenvalues, h, tolRate, s_0, t - t_now, N_SUBSTEPS_MIN, numericalErrorEstimate, &t_step, tmpKrylovVec1, &err_step);
			double s_0 = INITIAL_STEP_FRACTION * t_step;
			if(t_now == 0)
				findMaximalStepSize2(schurvector, eigenvalues, h, tolRate, s_0, t - t_now, N_SUBSTEPS_MIN, numericalErrorEstimate, &t_step, tmpKrylovVec1, &err_step, false);
			else
				findMaximalStepSize2(schurvector, eigenvalues, h, tolRate, s_0, t - t_now, N_SUBSTEPS_MIN, numericalErrorEstimate, &t_step, tmpKrylovVec1, &err_step, false);

            if (t_step <= 0)
            {
                std::cout << "Internal error: negative step size" << std::endl;
                exit(1);
            }
            cblas_zgemv(CblasColMajor, CblasNoTrans, Hsize, m, &one, V->values, Hsize, tmpKrylovVec1, 1, &zero, currentVec, 1);
        }
		//END STEP 2
        
		//STEP 3: sample observables in current time step, i.e. between 't_now' and 't_now + t_step'
		//Last condition of while-condition concerns last sampling point in case there was a numerical error in the addition of the time steps
        while (index_samples < n_samples && (t_sampling + samplingStep <= t_now + t_step || (index_samples == n_samples - 1 && std::abs(t_now + t_step - t) <= std::numeric_limits<double>::epsilon())))
        {
            t_sampling += samplingStep;
            //Determine state at t_sampling in the following
            cblas_zcopy(m, eigenvalues, 1, tmpKrylovVec1, 1);
            cblas_zdscal(m, (t_sampling - t_now), tmpKrylovVec1, 1);
            vzExp(m, tmpKrylovVec1,tmpKrylovVec2);
            //Now temp1 is no longer needed and can be reused
            cblas_zgemv(CblasColMajor, CblasConjTrans, m, m, &one, schurvector, m, e_1, 1, &zero, tmpKrylovVec1, 1);
            for(size_t i = 0; i != m; i++)
            {
            	tmpKrylovVec1[i] = tmpKrylovVec1[i]*tmpKrylovVec2[i];
            }
            //Now temp2 is no longer needed and can be reused
            cblas_zgemv(CblasColMajor, CblasNoTrans, m, m, &one, schurvector, m, tmpKrylovVec1, 1, &zero, tmpKrylovVec2, 1);
            cblas_zgemv(CblasColMajor, CblasNoTrans, Hsize, m, &one, V->values, Hsize, tmpKrylovVec2, 1, &zero, sampledState, 1);

            sample();
        }
        //END STEP 3
        
        t_now += t_step;
        err += err_step;
        
        //As a consistency check, determine if norm of state was preserved
        double nrm = cblas_dznrm2(Hsize, currentVec, 1);
        if (checkNorm)
        {
            if (fabs(nrm - 1) > tol)
            {
                std::cerr << "Error: computed norm is not inside specified tolerance." << std::endl;
                exit(1);
            }
        }
        
    }

    //Return result
    krylovReturn* ret = new krylovReturn(nbObservables, Hsize, n_samples);
    cblas_zcopy(Hsize, sampledState, 1, ret->evolvedState, 1);
    unsigned int nbResults;
    if (nbObservables != 0)
        nbResults = n_samples*nbObservables;
    else
        nbResults = n_samples *Hsize;
    cblas_zcopy(nbResults, samplings->values, 1, ret->sampling->values, 1);
    ret->n_steps = n_steps;
    ret->err = err;
    ret->dim = m;
    ret->nSamples = n_samples;
    
    delete[] eigenvalues;
    delete[] schurvector;
    delete V;
    delete H;
    destroyOptimizeInput();

    return ret;
}

/**
 * Perform the Arnoldi algorithm (simplified for an anti-Hermitian matrix) for the current state ('currentVec')
 * @param tolRate Maximal allowable error rate (for detection of happy breakdown)
 * @param HRet Returns the Hessenberg matrix
 * @param VRet Returns the corresponding transformation matrix
 * @param hRet Returns the element (m+1,m) of the Hessenberg matrix
 * @param mRet Returns the actual size of the Krylov subspace (important in case of happy breakdown)
 * @return false, is no happy breakdown has occured; true if happy breakdown has occured
 */
bool krylovTimeEvolver::arnoldiAlgorithm(double tolRate, matrix *HRet, matrix *VRet, double *hRet, size_t *mRet) {
	double normy = 0.;
	std::complex<double> negativeH;
	cblas_zcopy(Hsize, currentVec, 1, VRet->values, 1);

	for (size_t j = 0; j <= m - 1; j++) {
		sparse_status_t mklStatus = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, expFactor, *HamOpt,
				descriptor, (VRet->values) + j * Hsize, zero, tmpBlasVec);
        if(SPARSE_STATUS_SUCCESS != mklStatus)
        {
            std::cerr << "MKL error " << mklStatus << std::endl;
            exit(1);
        }
		if (j != 0) {
			negativeH = (-1.0) * *((HRet->values) + j - 1 + j * m);
			cblas_zaxpy(Hsize, &negativeH, VRet->values + (j - 1) * Hsize, 1,
					tmpBlasVec, 1);
		}
		cblas_zdotc_sub(Hsize, VRet->values + j * Hsize, 1, tmpBlasVec, 1,
				(HRet->values) + j + (j * m));
		negativeH = (-1.0) * *((HRet->values) + j + j * m);
		cblas_zaxpy(Hsize, &negativeH, VRet->values + j * Hsize, 1, tmpBlasVec,
				1);
		
		normy = cblas_dznrm2(Hsize, tmpBlasVec, 1);

		//Detection of lucky breakdown
		if (normy < tolRate) {
			std::cout << "***Lucky breakdown at Krylov dimension " << j + 1
					<< " *** " <<std::endl;
			*mRet = j + 1;
			*hRet = normy;
			return true;
		}
		//End detection of happy breakdown

		if (j + 1 != m) {
			HRet->values[j + (j + 1) * m].real(-normy);
			HRet->values[j + 1 + j * m].real(normy);
			std::complex<double> inverseNorm;
			inverseNorm.real(1.0 / normy);
			inverseNorm.imag(0.0);
			cblas_zscal(Hsize, &inverseNorm, tmpBlasVec, 1);
			cblas_zcopy(Hsize, tmpBlasVec, 1, (VRet->values) + Hsize * (j + 1),1);
		} else
			*hRet = normy;
	}
	*mRet = m;
	return false;
}

/*
* Compute the error integral to determine the analytic error of the krylov approximation
* @param a Start of integration
* @param b End of integration
* @param T Transformation matrix
* @param spectrumH Eigenvalue spectrum
* @param h Last entry of the Hessenberg matrix
*/
double krylovTimeEvolver::integrateError(double a, double b, std::complex<double>* T, std::complex<double>* spectrumH, double h)
{
	double error, L1;
	double termination = std::sqrt(std::numeric_limits<double>::epsilon());
	termination = 1.0/100.0;
	size_t level;
	auto f = [&](double x) {return h * std::abs(expKrylov(x, T, spectrumH)[m-1]); };
	double ret;
	size_t maxref = 8;



	//ret = integ.integrate(f, a, b, termination, &error, &L1, &level);

	ret = boost::math::quadrature::trapezoidal(f, a, b, termination, maxref, &error);

	//ret = boost::math::quadrature::gauss<double, 15>::integrate(f, a, b);

	//ret = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, a, b, 0, termination, &error);

	//boost::math::quadrature::tanh_sinh<double> intTanh(maxref);
	//ret = intTanh.integrate(f, a, b, termination, &error, &L1, &level);


	return ret;
}

/*
* Kernel of the error integral -- DELETE
*/
double krylovTimeEvolver::errorKernel(double t, std::complex<double>* T, std::complex<double>* spectrumH, double h)
{
	cblas_zcopy(m, spectrumH, 1, tmpintKernelExp1, 1);
	cblas_zdscal(m, t, tmpintKernelExp1, 1);
	vzExp(m, tmpintKernelExp1, tmpintKernelExp2);
	cblas_zgemv(CblasColMajor, CblasConjTrans, m, m, &one, T, m, e_1, 1, &zero, tmpintKernelT, 1);
	for (unsigned int i = 0; i != m; i++)
		tmpintKernelExp2[i] = tmpintKernelExp2[i] * tmpintKernelT[i];
	cblas_zgemv(CblasColMajor, CblasNoTrans, m, m, &one, T, m,
		tmpintKernelExp2, 1, &zero, tmpintKernelExp1, 1);

	double ret = h * std::abs(tmpintKernelExp1[m - 1]);
	return ret;
}

/*
* Calculate time evolution in Kyrlov space
*/
std::complex<double>* krylovTimeEvolver::expKrylov(double t, std::complex<double>* T, std::complex<double>* spectrumH)
{
	cblas_zcopy(m, spectrumH, 1, tmpintKernelExp1, 1);
	cblas_zdscal(m, t, tmpintKernelExp1, 1);
	vzExp(m, tmpintKernelExp1, tmpintKernelExp2);
	cblas_zgemv(CblasColMajor, CblasConjTrans, m, m, &one, T, m, e_1, 1, &zero, tmpintKernelT, 1);
	for (unsigned int i = 0; i != m; i++)
		tmpintKernelExp2[i] = tmpintKernelExp2[i] * tmpintKernelT[i];
	cblas_zgemv(CblasColMajor, CblasNoTrans, m, m, &one, T, m,
		tmpintKernelExp2, 1, &zero, tmpintKernelExp, 1);

	return tmpintKernelExp;

}