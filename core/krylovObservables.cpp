#include "krylovObservables.h"

using namespace TE;


/**
* Function to return the name of the observable
* @return Name of the observable
*/
std::string krylovBasicObservable::retName()
{
	return obs_name;
}

/**
* Function to return a pointer to the expectation value array
* @return Raw pointer to the array which stores the expectation values
*/
double* krylovBasicObservable::retExpectationValues()
{
    return expectationValues;
}

/**
* Function to return the number of samples taken
* @return Number of samples taken
*/
size_t krylovBasicObservable::retNumSamples()
{
    return numSamples;
}


/**
 * Initializes a raw array to requested size for storing the expectation values of the oberservable 
 * @param size number of entries in the result)
 */
void krylovBasicObservable::initializeResultArray(size_t size)
{	
	numSamples = size;
	expectationValues = new double[size];
}

/**
* Resets the index of sampling to zero. Note that the original values are not overwritten! 
* @return The size of the expectation value array, i.e. the number of expectation values that can be recorded.
*/
size_t krylovBasicObservable::resetResultArray()
{
    sampleIndex = 0;
    return numSamples;
}

/**
* Constructor for base class for observables, initializing the array for the expectation values with predefinied values
* @param name Name of the observable
* @param values Values to be stored in the expectation value array
*/
krylovBasicObservable::krylovBasicObservable(const std::string& name, std::vector<double> values): obs_name(name), dim(0), numSamples(values.size()), sampleIndex(values.size()), type(VOID_TYPE_OBS), expectationValues(nullptr)
{
    expectationValues = new double[numSamples];

    for (int i = 0; i != numSamples; i++) {
        expectationValues[i] = values[i];
    }
}

/**
* Default deconstructor
*/
krylovBasicObservable::~krylovBasicObservable()
{
	if (numSamples > 0)
		delete[] expectationValues;
}


/**
* Function to return the type of the observable. 
* @return Type of the observable
*/
obsType krylovBasicObservable::retType()
{
	return type;
}


/**
* Constructor for (dense) matrix observables 
* @param name Name of the observable
* @param obser (dense) Matrix representation of the observable 
*/
krylovMatrixObservable::krylovMatrixObservable(const std::string& name, std::unique_ptr<matrix> obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = MATRIX_TYPE_OBS;
	if (dim > 0)
		tmpBlasVec = new std::complex<double>[dim]; //array for storing temporary intermediate values

	obs = std::move(obser);

}

/**
* Destructor for (dense) matrix observables
*/
krylovMatrixObservable::~krylovMatrixObservable()
{
	if (dim > 0)
		delete[] tmpBlasVec;
}


/**
 * Computes expectation value of a dense matrix observable for a given quantum state
 * @param vec Quantum state vector
 * @param len Length of state vector
 */
std::complex<double> krylovMatrixObservable::expectation(std::complex<double>* vec, int len) //requires testing
{
    if (sampleIndex >= numSamples)
    {
        std::cerr << "Too many samples." << std::endl;
        exit(1);
    }
	if (len != dim)
	{
		std::cerr << "Incompatible dimensions" << std::endl;
		exit(1);
	}
	std::complex<double> observall;
	cblas_zgemv(CblasColMajor, CblasNoTrans, dim, dim, &one, obs->values, dim, vec, 1, &zero, tmpBlasVec, 1);
	cblas_zdotc_sub(len, vec, 1, tmpBlasVec, 1, &observall);

	expectationValues[sampleIndex] = observall.real();
	sampleIndex++;
	return observall;
}

/**
* Constructor for (sparse) matrix observables
* @param name Name of the observable
* @param obser (sparse) Matrix representation of the observable
*/
krylovSpMatrixObservable::krylovSpMatrixObservable(const std::string& name, std::unique_ptr<smatrix> obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = SPARSE_MATRIX_TYPE_OBS;
	obs = std::move(obser); 
	obs->initialize();

	tmpBlasVec = new std::complex<double>[dim];
}

/**
* Destructor for (sparse) matrix observables
*/
krylovSpMatrixObservable::~krylovSpMatrixObservable()
{
	if (dim > 0)
	{
		delete[] tmpBlasVec;
	}
}


/**
 * Computes expectation value of a sparse matrix observable for a given quantum state
 * @param vec Quantum state vector
 * @param len Length of state vector
 * @return computed expectation value <Obs|state|Obs>. Since we do technically not require the observable to be Hermitian the value might be complex. Note that the stored value only contains the real part. 
 */
std::complex<double> krylovSpMatrixObservable::expectation(std::complex<double>* vec, int len)
{
    if (sampleIndex >= numSamples)
    {
        std::cerr << "Too many samples." << std::endl;
        exit(1);
    }
	if (len != dim)
	{
		std::cerr << "Incompatible dimensions" << std::endl;
		exit(1);
	}

	std::complex<double> observall;
	obs->spMV(one, vec, tmpBlasVec);
	cblas_zdotc_sub(len, vec, 1, tmpBlasVec, 1, &observall);
	expectationValues[sampleIndex] = observall.real();
	sampleIndex++;

	return observall;

}

/**
* Constructor for vector observables
* @param name Name of the observable
* @param obser Pointer to an complex array representing the vector
* @param len Length of the array 
*/
krylovVectorObservable::krylovVectorObservable(const std::string& name, std::complex<double>* obser, size_t len) : krylovBasicObservable(name)
{
	dim = len;
	type = VECTOR_TYPE_OBS;
	obs = std::make_unique<std::complex<double>[]>(len);
	cblas_zcopy(dim, obser, 1, obs.get(), 1);
}


/**
 * Computes expectation value of a vector observable for a given quantum state or equivalently the absolute square of a projection of a complexvector onto the state vector
 * @param vec Quantum state vector
 * @param len Length of state vector
 */
std::complex<double> krylovVectorObservable::expectation(std::complex<double>* vec, int len)
{
    if (sampleIndex >= numSamples)
    {
        std::cerr << "Too many samples." << std::endl;
        exit(1);
    }
	if (len != dim)
	{
		std::cerr << "Incompatible dimensions" << std::endl;
		exit(1);
	}
	std::complex<double> observall;
	cblas_zdotc_sub(len, vec, 1, obs.get(), 1, &observall);
	std::complex<double> observallreturn;
	observallreturn.imag(0);
	observallreturn.real(std::norm(observall));
	expectationValues[sampleIndex] = observallreturn.real();
	sampleIndex++;
	return observallreturn; //returns the squared magnitued
}

/**
* Constructor for the output observable class. Note that this class is only for file output. 
* @param name Name of the observable
* @param values Vector of values that should be stored in the output observable
*/
krylovOutputObservable::krylovOutputObservable(const std::string& name, std::vector<double> values) : krylovBasicObservable(name)
{
    this->numSamples = values.size();
    initializeResultArray(numSamples);
    for (int i = 0; i != numSamples; i++) {
        expectationValues[i] = values[i];
    }
}

/**
* Throw runtime error when output observable is tried to be sampled. 
*/
std::complex<double> krylovOutputObservable::expectation(std::complex<double>* vec, int len)
{
    throw std::runtime_error("krylovOutputObservable should not call `expectation` function.");
    return 0;
}


/**
* Default exception message
* @return char array of the message
*/
const char* requestStopException::what() const throw()
{
    return "Observable requested termination of time evolution.";
}


