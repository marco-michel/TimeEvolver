#include "krylovObservables.h"






using namespace TE;


std::string krylovBasicObservable::retName()
{
	return obs_name;
}

double* krylovBasicObservable::retexpectationValues()
{
    return expectationValues;
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
 * Writes sampled expectation values of a list of observables to file
 * @param obs_list List of observables with sampled expectation values
 * @param name Requested filename
 */
void krylovBasicObservable::saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>> &obs_list,  parameter_list &para, const std::string &name)
{
    std::string outputFileName = name;

    parameter_list::iterator paraIter;
    std::vector<std::unique_ptr<krylovBasicObservable>>::const_iterator obsIter;

    //Create string with parameter and its value
    for (paraIter = para.begin(); paraIter != para.end(); paraIter++)
        outputFileName += "_" + (*paraIter)->getName() + (*paraIter)->getData();

    obsIter = obs_list.begin();
    size_t nbOutputParameters = para.size();


#ifdef USE_HDF
    outputFileName += ".h5";
    H5File fileHh(outputFileName.c_str(), H5F_ACC_TRUNC);
    DataSet dataset;

    for (unsigned int j = 0; j < obs_list.size() && obsIter != obs_list.end(); j++)
    {
        int NX = (int)(*obsIter)->numSamples;
        const int RANK = 1;
        hsize_t dimsf[RANK];
        dimsf[0] = NX;
        DataSpace dataspace(RANK, dimsf);
        FloatType datatype(PredType::NATIVE_DOUBLE);
        datatype.setOrder(H5T_ORDER_LE);
        std::string observableName = (*obsIter)->obs_name;
        dataset = fileHh.createDataSet(observableName.c_str(), datatype, dataspace);
        dataset.write((*obsIter)->expectationValues, PredType::NATIVE_DOUBLE);

        //write  parameters as attributes to file
        hsize_t dims[1] = { 1 };

        DataSpace** attr_dataspace = new DataSpace * [nbOutputParameters];
        Attribute** attributes = new Attribute * [nbOutputParameters];

        for (unsigned int i = 0; i < nbOutputParameters; ++i)
        {
            attr_dataspace[i] = new DataSpace(1, dims);
            attributes[i] = new Attribute();
        }

        int counter = 0;
        for (paraIter = para.begin(); paraIter != para.end(); paraIter++, counter++)
        {
            if ((*paraIter)->isDouble())
            {
                double paraValue = dynamic_cast<typedParameter<double>&>(*(*paraIter)).getValue();
                *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_DOUBLE, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_DOUBLE, &paraValue);
            }
            else if ((*paraIter)->isInt())
            {
                int paraValue = dynamic_cast<typedParameter<int>&>(*(*paraIter)).getValue();
                *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_INT, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_INT, &paraValue);
            }
            else if ((*paraIter)->isBool())
            {
                int paraValue = dynamic_cast<typedParameter<bool>&>(*(*paraIter)).getValue();
                *attributes[counter] = dataset.createAttribute((*paraIter)->getName(), PredType::NATIVE_INT, *attr_dataspace[counter]); attributes[counter]->write(PredType::NATIVE_INT, &paraValue);
            }
        }

        for (unsigned int i = 0; i < nbOutputParameters; ++i)
        {
            delete attr_dataspace[i];
            delete attributes[i];
        }

        delete[] attr_dataspace;
        delete[] attributes;

        dataset.close();
        obsIter++;

    }
    //If not write data to simple csv files
#else

for (; obsIter != obs_list.end(); obsIter++)
{
    std::string fileNameCSV = outputFileName + (*obsIter)->retName() + ".csv";
    std::ofstream outputfile;
    outputfile.open(fileNameCSV);

    
    for (int i = 0; i != (*obsIter)->numSamples - 1; i++)
        outputfile << (*obsIter)->expectationValues[i] << ", ";
    outputfile << (*obsIter)->expectationValues[((*obsIter)->numSamples) - 1];
    outputfile.close();
}

#endif

}



krylovBasicObservable::~krylovBasicObservable()
{
	if (numSamples > 0)
		delete[] expectationValues;
}

obsType krylovBasicObservable::retType()
{
	return type;
}


krylovMatrixObservable::krylovMatrixObservable(const std::string& name, std::unique_ptr<matrix> obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = MATRIX_TYPE_OBS;
	if (dim > 0)
		tmpBlasVec = new std::complex<double>[dim]; //array for storing temporary intermediate values

	obs = std::move(obser);

}

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

krylovSpMatrixObservable::krylovSpMatrixObservable(const std::string& name, std::unique_ptr<smatrix> obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = SPARSE_MATRIX_TYPE_OBS;
	obs = std::move(obser); 
	obs->initialize();

	tmpBlasVec = new std::complex<double>[dim];
}

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


krylovVectorObservable::krylovVectorObservable(const std::string& name, std::complex<double>* obser, size_t len) : krylovBasicObservable(name)
{
	dim = len;
	type = VECTOR_TYPE_OBS;
	obs = std::make_unique<std::complex<double>[]>(len);
	cblas_zcopy(dim, obser, 1, obs.get(), 1);
}


/**
 * Computes expectation value of a vector observable for a given quantum state or equivalently a projection of a complexvector onto a state vector
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
