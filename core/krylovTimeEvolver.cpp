/** Brief summary :
 For given Hamiltonian, time evolve a given initial state, while also saving results at intermediate time steps.
 This is done in such a way that at each step, the norm difference of the state obtained by numerical time evolution deviates from the true result at most by a prescribed error bound.
 The numerical routine is based on iteratively finding a Krylov-subspaces, on which time evolution can be perform for a small time step.
 The time-evolved vector then becomes the initial vector for constructing the subsequent Krylov-subspace.
 */

#include "krylovTimeEvolver.h"

using namespace TE;


/**
 * All data required for the numerical time evolution is needed.
 * @param t The time interval over which the state should be time evolved.
 * @param v The initial and normalized state that should be time evolved
 * @param samplingStep The time interval after which the values of the observables should be determined
 * @param tol The maximal admissible error (norm difference between result of numerical and true time evolution)
 * @param m The size of the Krylov subspaces
 * @param observables The vector of observables that are to be sampled
 * @param Ham The full Hamiltonian
 * @param expFactor The scalar factor multiplying the Hamiltonian in the time evolution exp(-i expFactor H t)
 * @param fastIntegration Whether a faster but less accurate method for evaluating the error integral should be used (default value: false)
 * @param progressBar Whether or not to show a progressbar in the terminal
 */
krylovTimeEvolver::krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, double tol, int mm, std::vector<std::unique_ptr<krylovBasicObservable>>  observables, std::unique_ptr<smatrix> HamIn, double expFactor, bool fastIntegration, bool progressBar)
{


	this->t = t;  this->samplingStep = samplingStep; this->tol = tol;  this->progressBar = progressBar;
	this->Ham = std::move(HamIn); this->expFactor = std::complex<double>(0,-1*expFactor); this->fastIntegration = fastIntegration; this->nbObservables = (int)observables.size();

	matrixNorm = Ham->norm1();
	this->Hsize = Ham->m;
	this->m = std::min<size_t>(mm, Hsize);

	stop_printing = false;

	if (Hsize == 0)
	{
		logger.log_message(krylovLogger::FATAL, "Invalid Hilbertspace dimension");
		exit(1);
	}

	if (std::abs(cblas_dznrm2(Hsize, v, 1) - 1.0) > tol) {
		logger.log_message(krylovLogger::FATAL, "Initial vector is not normalized");
		exit(1);
	}

	//Perform initializing steps for the matrix representation tailored to the installed (s)BLAS library
	Ham->initialize();

	//Numerical integration terminates if error*L1 < termination
	termination = 1e-3;

	//number of sampling steps
	n_samples = (size_t)floor(t / samplingStep) + 1;

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
	tmpintKernelExp3 = new std::complex<double>[m];
	tmpintKernelT = new std::complex<double>[m];

	index_samples = 0;
	e_1 = new std::complex<double>[m];
	e_1[0].real(1);

	obsVector = std::move(observables);
	for (auto iter = obsVector.begin(); iter != obsVector.end(); iter++)
		(*iter)->initializeResultArray(n_samples);

}


krylovTimeEvolver::krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, std::vector<std::unique_ptr<krylovBasicObservable>> observables, std::unique_ptr<smatrix> Ham) : krylovTimeEvolver(t, v, samplingStep, 1e-6, 40, std::move(observables), std::move(Ham),
	1.0, false, false){}

/**
* Destructor
*/
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
	delete[] tmpintKernelExp3;
	delete[] tmpintKernelT;
    delete[] e_1;
	obsVector.clear();
}



/**
 * Computes and saves values of observables for current state ('sampledState')
 */
void krylovTimeEvolver::sample() {
	int i = 0;
	for (auto obsIter = obsVector.begin(); obsIter != obsVector.end(); obsIter++, i++)
	{
		(*obsIter)->expectation(sampledState, Hsize);
	}
	index_samples++;
}

/**
* Sets the severity level for logging
* @param level: Specific level, e.g. DEBUG, WARNING, ERROR, ... 
*/
void krylovTimeEvolver::changeLogLevel(krylovLogger::loggingLevel level)
{
	logger.set_loggingLevel(level);
}

/**
 * Given a Krylov subspace, this function finds out how far (i.e. for what time step) one can use it without exceeding the prescribed error bound. To this end, it computes an error bound for different possible time steps. This function only operates in the Krylov subspace, i.e. with vectors and matrices of dimension m
 * @param T The orthogonal transformation matrix for the Hamiltonian in the Krylov subspace
 * @param spectrumH The eigenvalues of the Hamiltonian in the Krylov subspace
 * @param h The last element of the Arnoldi algorithm (it's needed for the computation of the error)
 * @param tolRate The maximal admissible error rate (i.e. error per time)
 * @param t_stepSuggestion The initial substep size
 * @param t_step_max The maximal stepsize (in case tolRate is not exceeded)
 * @param n_substeps The number of substeps (this determines how finely the program tries to increase the time step)
 * @param numericalErrorEstimate The numerical error associated to the Krylov subspace (this sets a lower bound for the error)
 * @param increaseStep Whether or not the function should first try to increase the step size (useful if no good suggestion was provided)
 * @param t_stepRet Returns the time step
 * @param w_KrylovRet Returns the Krylov vector after the time step
 * @param err_stepRet Returns the error of the time step
 * @return Return errorcode (0=Success, 1=estimate of roundoff errors are larger than analytical error, 20=error of integration larger than requested termination tolerance, multiple errors are indicated as the sum of respective error codes)
 */
int krylovTimeEvolver::findMaximalStepSize(std::complex<double>* T, std::complex<double>* spectrumH, double h, double tolRate, double t_stepSuggestion, double t_step_max, int n_substeps, double numericalErrorEstimate, bool increaseStep, double* t_stepRet, std::complex<double>* w_KrylovRet, double* err_stepRet)
{
	//Maximal number of substepreductions to meet tolerance
	unsigned int GO_MAX = 20;
	unsigned int nbReductions = 0;

	double t_step = t_stepSuggestion;
	int statusCode = 0;
	bool numericalIntegrationSuccessful = true; //Did numerical integration converge to adequate tolerance?
	bool skipSubsteps = false; //needed if t_step > t_step_max


	if (t_step > t_step_max)
	{
		t_step = t_step_max; 
		skipSubsteps = true; //there is no need to increase the step size anylonger 
	}


	//Integrate to get a first error estimate 
	double err_step = integrateError(0, t_step, T, spectrumH, h, integrationMethodLong, tolRate, numericalIntegrationSuccessful);
	
	//Increasing step size. Used in the first go through when there is no good guess for the step size. 
	if (increaseStep && !skipSubsteps)
	{
		double err_step_new = err_step;
		double t_step_new = t_step;
		while (err_step_new < tolRate * t_step_new && t_step < t_step_max)
		{
			t_step = t_step_new;
			err_step = err_step_new;
			if(t_step>= t_step_max){
				skipSubsteps = true;
				break;
			}
			t_step_new *= 2;
			if(t_step_new > t_step_max)
			{
				t_step_new = t_step_max;
			}
			err_step_new = integrateError(0, t_step_new, T, spectrumH, h, integrationMethodLong, tolRate, numericalIntegrationSuccessful);	
		}
	} 
	
	//Reduce step_size in case the error is not within requested tolerance
	while (err_step > tolRate * t_step)
	{
		nbReductions++;
		t_step = t_step / 2.0;
		err_step = integrateError(0, t_step, T, spectrumH, h, integrationMethodLong, tolRate, numericalIntegrationSuccessful);
		if (nbReductions == GO_MAX)
		{
			logger.log_message(krylovLogger::FATAL, "Error: No small enough time step found to meet tolerance requirements.");
			exit(1);
		}
		skipSubsteps = false;
	}
	
	double s = t_step / n_substeps;
	int n_s = 0;
	double deltaError = 0;
	
	//Increase step size in small substeps to use the Krylov space as long as possible
	if (!skipSubsteps)
	{
		while (err_step + deltaError < tolRate * (t_step + n_s * s))
		{
			err_step += deltaError;
			n_s++;
			if(t_step +(n_s-1)*s >= t_step_max)
			{
				break;
			} 
			else
			{
				deltaError = integrateError(t_step + (n_s-1) * s, t_step + n_s * s, T, spectrumH, h, integrationMethodShort, tolRate, numericalIntegrationSuccessful);
			}
		}

		t_step += (n_s - 1) * s;
	}

	if (err_step < numericalErrorEstimate)  //Estimate of numerical error is larger than analytic error
	{
		if(!skipSubsteps || err_step > 0.1 * tolRate * t_step) //Warning can only be skipped if it's the last Krylov space and the error budget was barely exploited. 
			statusCode += 1;
	}
		

	if (!numericalIntegrationSuccessful)
		statusCode += 20; 	//Numerical integration can not be trusted, nullifing the error bound

	*t_stepRet = t_step;
	*err_stepRet = err_step;

	//Compute state in Krylov space at t=t_step and return it
	std::complex<double>* w_Krylov = expKrylov(t_step, T, spectrumH);
	cblas_zcopy(m, w_Krylov, 1, w_KrylovRet, 1);

	return statusCode;

}


/**
 * Main rountine: performs time evolution according to the data provided in the constructor
 * @return Returns in particular the samplings for each observable, the state after time evolution and an upper bound on the error
 */
krylovReturn* krylovTimeEvolver::timeEvolve()
{
	//Constants
	//Number of substeps per time step (determines how finely the program tries to increase the maximal timestep per Krylov space)
	int N_SUBSTEPS = 50;
	//After each time step, the optimal step size is computed. To avoid substep reduction because of a too small number of substeps, the optimal step size is multiplied by this number
	double INITIAL_STEP_FRACTION = 0.97;

	if (fastIntegration)//Use Gauss integration
	{
		integrationMethodLong = 0;
		integrationMethodShort = 1;
	}
	else//Use double exponential integration
	{
		integrationMethodLong = integrationMethodShort = 2;
	}

	if (fastIntegration)
		logger.log_message(krylovLogger::WARNING, "Please note that a less accurate method for evaluating the error integral was used. We recommend recomputing with accurate integration to establish the validity of the result.");



	//Estimate for initial step size
	double t_step = 1.0 / 10.0;
	//Admissible error rate
	double tolRate = tol / t;
	//Estimate of error due to limited precision of numerical operations per Kyrlov space
	double numericalErrorEstimate = Hsize * matrixNorm * std::numeric_limits<double>::epsilon();
	//Total esimate of round-off errors
	double numericalErrorEstimateTotal = 0;

	//Flag indicating if a lucky breakdown has occured
	bool dummy_hbd = false;
	//In case of lucky breakdown, size of Krylov space
	size_t m_hbd;
	//Monitoring if something went wrong in findSubstep
	int errorCodeFindSubstep = 0;
	//Count number of Kyrlov spaces with error bound failures
	int nbErrRoundOff = 0;
	int nbErrInt = 0;
	//Were there multible errors thrown
	int nbErrors = 0;

	//Record observables for initial state
	try {
		sample();
	}
	catch (requestStopException& e)
	{
		return generateReturn();
	}


	//Hessenberg matrix
	matrix* H = new matrix(m, m);
	//Corresponding transformation matrix
	matrix* V = new matrix(Hsize, m);
	//The (m+1,m) element of Hessenberg matrix (needed for computation of error)
	double h = 0;
	//Eigenvalues of Hessenberg matrix
	std::complex<double>* eigenvalues = new std::complex<double>[m];
	//Eigenvectors of Hessenberg matrix
	std::complex<double>* schurvector = new std::complex<double>[m * m];

	//Init observable vector before first sample
	



	//Start progressBar thread
	if (progressBar == true)
			pBThread = std::thread(&krylovTimeEvolver::progressBarThread, this);
	
    //Main loop
    while (index_samples < n_samples)
    {
		n_steps++;
		//t_step = fmin(t - t_now, t_step);
		double err_step = 0;

		//STEP 1: Construct Krylov subspace using Arnoldi algorithm
		dummy_hbd = arnoldiAlgorithm(tolRate, H, V, &h, &m_hbd);

		//Some special adjustments in case of a lucky breakdown, i.e. when projection in Krylov-subspace of dimension m_hbd <= m is exact (within numerical uncertainty)
		//In particular, the time step of the current Krylov space can be arbitarily large in this case
		if (dummy_hbd) 
		{
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

			logger.log_message(krylovLogger::INFO, "***Lucky breakdown at Krylov dimension" + std::to_string(m) + " *** ");

			statusCode = 1;
		}
		//Finally diagonalize Hessenberg matrix H (since it will be exponentiated many times)
		int infocheck = TE_zhseqr(m, H->values, eigenvalues, schurvector);
		if (infocheck != 0) 
		{
			logger.log_message(krylovLogger::FATAL, "Internal error: LAPACK error " + infocheck);
			exit(1);
		}
		//END STEP 1
        
        //STEP 2: find maximal time step ('t_step') for which current Krylov subspace can be used without violating error bound
		if (dummy_hbd == false)
        {
			double s_0 = INITIAL_STEP_FRACTION * t_step;
			if(t_now == 0)
				errorCodeFindSubstep = findMaximalStepSize(schurvector, eigenvalues, h, tolRate, s_0, t - t_now, N_SUBSTEPS, numericalErrorEstimate, true, &t_step, tmpKrylovVec1, &err_step);
			else
				errorCodeFindSubstep = findMaximalStepSize(schurvector, eigenvalues, h, tolRate, s_0, t - t_now, N_SUBSTEPS, numericalErrorEstimate, false, &t_step, tmpKrylovVec1, &err_step);

            cblas_zgemv(CblasColMajor, CblasNoTrans, Hsize, m, &one, V->values, Hsize, tmpKrylovVec1, 1, &zero, currentVec, 1);
        }

		//Error handling of findMaximalStepSize
		switch (errorCodeFindSubstep)
		{
		case 0:
			break; //There were no problems detected
		case 1:
			nbErrRoundOff++; //Numerical error was larger than analytic error
			break;
		case 20:
			nbErrInt++; //Requested accuracy of numerical integration could not be met  
			break;
		case 21:
			nbErrInt++; nbErrRoundOff++;
			break;
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
			expV(m, tmpKrylovVec1, tmpKrylovVec2);

            //Now temp1 is no longer needed and can be reused
            cblas_zgemv(CblasColMajor, CblasConjTrans, m, m, &one, schurvector, m, e_1, 1, &zero, tmpKrylovVec1, 1);
            for(size_t i = 0; i != m; i++)
            {
            	tmpKrylovVec1[i] = tmpKrylovVec1[i]*tmpKrylovVec2[i];
            }
            //Now temp2 is no longer needed and can be reused
            cblas_zgemv(CblasColMajor, CblasNoTrans, m, m, &one, schurvector, m, tmpKrylovVec1, 1, &zero, tmpKrylovVec2, 1);
            cblas_zgemv(CblasColMajor, CblasNoTrans, Hsize, m, &one, V->values, Hsize, tmpKrylovVec2, 1, &zero, sampledState, 1);

			try {
				sample();
			} catch (requestStopException& e)
			{
				delete[] eigenvalues; delete[] schurvector;
				delete V; delete H;
				if (progressBar){
					stop_printing = true;
					pBThread.join();
				}
				statusCode = 1000;
				return generateReturn();
			}
        }
        //END STEP 3
        
        t_now += t_step;
        err += err_step;
    }


	//As a consistency check, determine if norm of state was preserved
	if (std::abs(cblas_dznrm2(Hsize, sampledState, 1) - 1.0) > tol)
	{
		logger.log_message(krylovLogger::WARNING, "CRITICAL WARNING: Norm of state vector is not inside specified tolerance. \n \
				THE DESIRED ERROR BOUND WILL LIKELY BE VIOLATED.");
		nbErrors++;
		statusCode = 30;
	}


	if (progressBar)
		pBThread.join();

	//Output of potential errors
	numericalErrorEstimateTotal = numericalErrorEstimate * n_steps; //Rough estimate of total round-of error
	if (nbErrRoundOff != 0)
	{
		logger.log_message(krylovLogger::WARNING, std::string("CRITICAL WARNING: The computed error bound was smaller than the estimate of the numerical error in ") + (std::to_string(nbErrRoundOff) + " Krylov spaces.) \n \
			THE DESIRED ERROR BOUND WILL LIKELY BE VIOLATED. \n \
			Restart with bigger error bound or smaller time."));
		logger.log_message(krylovLogger::WARNING, "The total computed analytic error is: " + (std::stringstream() << err).str());
		logger.log_message(krylovLogger::WARNING, "The total numerical error estimate for all Krylov spaces is " + (std::stringstream() << numericalErrorEstimateTotal).str());

		nbErrors++;
		statusCode = 10;
	}
	if (nbErrInt != 0)
	{
		logger.log_message(krylovLogger::WARNING, "CRITICAL WARNING: The error of numerical integration did not meet the required accuarcy in " + std::to_string(nbErrInt) + " Krylov spaces. \
			THE DESIRED ERROR BOUND MAY BE VIOLATED.");

		nbErrors++;
		statusCode = 20;
	}
	if (numericalErrorEstimateTotal > err && nbErrRoundOff == 0) //No need to output warning twice. Only relevant for very simple systems for which 1 Krylov space is sufficient /  lucky breakdown
	{
		if (numericalErrorEstimateTotal > tol)
		{
			logger.log_message(krylovLogger::WARNING, "CRITICAL WARNING: The numerical error estimate " + (std::stringstream() << numericalErrorEstimateTotal).str() + " was larger than the requested tolerance " + (std::stringstream() << tol).str() + " \
				THE DESIRED ERROR BOUND WILL LIKELY BE VIOLATED.");

			statusCode = 11;
			nbErrors++;
		}
		else
		{
			
			logger.log_message(krylovLogger::WARNING, "Warning: Analytic error " + (std::stringstream() << err).str() + " was smaller than the estimate of the numerical error " + (std::stringstream() << numericalErrorEstimateTotal).str() + ". \
				The computed error is not accurate.");
			
			if (statusCode != 1) //Don't overwrite Lucky breakdown status
				statusCode = 2;
		}
	}

	if (nbErrors > 1)
	{
		logger.log_message(krylovLogger::WARNING, "There were multiple critical warnings raised during evaluation. Please see ouput above for further details.");
		statusCode = 100;
	}

    delete[] eigenvalues;
    delete[] schurvector;
    delete V;
    delete H;

    return generateReturn();
}

/**
 * Perform the Arnoldi algorithm (simplified for an anti-Hermitian matrix) for the current state ('currentVec')
 * @param tolRate Maximal allowable error rate (for detection of lucky breakdown)
 * @param HRet Returns the Hessenberg matrix
 * @param VRet Returns the corresponding transformation matrix
 * @param hRet Returns the element (m+1,m) of the Hessenberg matrix
 * @param mRet Returns the actual size of the Krylov subspace (important in case of lucky breakdown)
 * @return false, is no lucky breakdown has occured; true if lucky breakdown has occured
 */
bool krylovTimeEvolver::arnoldiAlgorithm(double tolRate, TE::matrix *HRet, TE::matrix *VRet, double *hRet, size_t *mRet) {
	double normy = 0.;
	std::complex<double> negativeH;
	cblas_zcopy(Hsize, currentVec, 1, VRet->values, 1);

	for (size_t j = 0; j <= m - 1; j++) {

		int spStatus = Ham->spMV(expFactor, (VRet->values) + j * Hsize, tmpBlasVec);

		//sparse_status_t mklStatus2 = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, expFactor, *HamOpt, descriptor, (VRet->values) + j * Hsize, zero, tmpBlasVec);

        if(spStatus != 0)
        {
			logger.log_message(krylovLogger::FATAL, "spMV error ");
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
			*mRet = j + 1;
			*hRet = normy;
			return true;
		}
		//End detection of lucky breakdown

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

/**
* Compute the error integral to determine the analytic error of the krylov approximation
* @param a Start of integration
* @param b End of integration
* @param T Transformation matrix
* @param spectrumH Eigenvalue spectrum
* @param h Last entry of the Hessenberg matrix
* @param method 0 stands for Gauss with 15 abscissa, 1 stands for Gauss with 7 abscissa, 2 stands for an adaptive sinh-tanh method
* @param tolRate maximal admissible error per time (used to compare scale of integral against it; if integral is very small, accuracy of integration will not be monitored)
* @param successful logical and-conjuction of input value and bool indicating whether numerical integration converged to sufficient accuracy (accuracy is not monitored in case of Gauss-integration so in this case the input value is always returned)
*/
double krylovTimeEvolver::integrateError(double a, double b, std::complex<double>* T, std::complex<double>* spectrumH, double h, int method, double tolRate, bool& successful)
{
	double error, L1;
	double ret;
	bool success = true;

	//Define Integrand as a lambda function
	auto f = [&](double x) {return h * std::abs(expKrylov(x, T, spectrumH)[m - 1]); };

	if(method == 0)
		ret = boost::math::quadrature::gauss<double, 15>::integrate(f, a, b); //Gauss-Legendre quadrature with 15 abscissa 
	else if(method == 1)
		ret = boost::math::quadrature::gauss<double, 7>::integrate(f, a, b); //Gauss-Legendre quadrature with 7 abscissa 
	else if (method == 2)
	{
		try {
			ret = integ.integrate(f, a, b, termination, &error, &L1); //Double exponential integration
		}
		catch (std::exception const& x) { //boost integration routine throw exception
			logger.log_message(krylovLogger::WARNING, boost::diagnostic_information(x));
			ret = this->tol * 10; //ignore the error in the integration by setting error to fail the following check to invoke substep reduction
		}
		//Check if requested accuracy was achieved; no check is needed if the contribution of the integral computed here to the absolute error integral is small
		if (error / ret > termination && ret > termination *(b-a) * tolRate)
			success = false;
	}
	else
	{
		logger.log_message(krylovLogger::FATAL, "Internal error: No method of integration selected");
		exit(1);
	}
	
	successful = successful && success;
	return ret;
}


/**
* Print progress status in terminal
* @param prog fraction of completion
*/
void krylovTimeEvolver::printProgress(float prog)
{
	std::cout << "[";
	int pos = (int) (pBarWidth * prog);
	for (int i = 0; i < pBarWidth; ++i)
	{
		if (i <= pos) std::cout << "|";
		else std::cout << " ";
	}
	std::cout << "] " << int(prog * 100.0) << " %\r";
	std::cout.flush();
}


/**
* ProgressBar thread
*/
void krylovTimeEvolver::progressBarThread()
{
	float prog;
	while (!stop_printing) {
		prog = static_cast<float>(index_samples) / n_samples;
		printProgress(prog);
		if (prog == 1) {
			printProgress(prog);
			break;
		}
	std::this_thread::sleep_for(std::chrono::milliseconds(500));
	}
	return;
}


/**
* Initialize member variables of return structure
*/
krylovReturn* krylovTimeEvolver::generateReturn()
{
	krylovReturn* ret = new krylovReturn(Hsize, statusCode);
	cblas_zcopy(Hsize, sampledState, 1, ret->evolvedState, 1);
	ret->n_steps = n_steps;
	ret->err = err;
	ret->dim = Hsize;
	ret->krylovDim = m;
	ret->observableList = std::move(obsVector);
	ret->hamiltonianMatrix = std::move(Ham);
	ret->evolvedTime = t_now;
	ret->numSamples = index_samples;

	return ret;
}


/**
* Calculate time evolution in Kyrlov space
* @param t Time
* @param T Transformation matrix
* @param spectrumH Eigenvalues
*/
std::complex<double>* krylovTimeEvolver::expKrylov(double t, std::complex<double>* T, std::complex<double>* spectrumH)
{
	cblas_zcopy(m, spectrumH, 1, tmpintKernelExp1, 1);
	
	//Exponenting scaled eigenvalues
	cblas_zdscal(m, t, tmpintKernelExp1, 1);
	expV(m, tmpintKernelExp1, tmpintKernelExp2);

	
	//Rotate basis back to Kyrlovspace
	cblas_zgemv(CblasColMajor, CblasConjTrans, m, m, &one, T, m, e_1, 1, &zero, tmpintKernelT, 1);
	mulV(m, tmpintKernelExp2, tmpintKernelT, tmpintKernelExp3);
	cblas_zgemv(CblasColMajor, CblasNoTrans, m, m, &one, T, m,
		tmpintKernelExp3, 1, &zero, tmpintKernelExp, 1);

	return tmpintKernelExp;

}
