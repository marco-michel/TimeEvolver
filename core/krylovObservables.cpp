
#include <iostream>
#include <complex>
#include <memory>

#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t
#include <mkl.h>
#elif defined USE_OPENBLAS
#include <cblas.h>
#endif

#include "matrixDataTypes.h"
#include "krylovObservables.h"
#include "krylovHelper.h"


using namespace TE;


std::string krylovBasicObservable::retName()
{
	return obs_name;
}

void krylovBasicObservable::initializeResultArray(size_t size)
{	
	numSamples = size;
	expectationValues = new double[size];
}



void krylovBasicObservable::saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>> &obs_list,  parameter_list &para, const std::string &name)
{
    std::string outputFileName = name;

    parameter_list::iterator paraIter;
    std::vector<std::unique_ptr<krylovBasicObservable>>::const_iterator obsIter;

    //Create string with parameter and its value
    for (paraIter = para.begin(); paraIter != para.end(); paraIter++)
        outputFileName += "_" + (*paraIter)->getName() + (*paraIter)->getData();

#ifdef USE_HDF

    outputFileName += ".h5";
    H5File fileHh(outputFileName.c_str(), H5F_ACC_TRUNC);
    DataSet dataset;

    size_t nbOutputParameters = para.size();

    obsIter = obs_list.begin();

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

for (int j = 0; j != nbObservables && obsIter != obs_list.end(); j++, obsIter++)
{
    std::string fileNameCSV = outputFileName + (*obsIter)->getName() + ".csv";
    std::ofstream outputfile;
    outputfile.open(fileNameCSV);
    for (int i = 0; i != (results->nSamples) - 1; i++)
        outputfile << nicelySorted[j][i] << ", ";
    outputfile << nicelySorted[j][(results->nSamples) - 1];
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


krylovMatrixObservable::krylovMatrixObservable(const std::string& name, matrix* obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = MATRIX_TYPE_OBS;
	if (dim > 0)
		tmpBlasVec = new std::complex<double>[dim];

	obs = std::make_unique<matrix>(*obser);

}

krylovMatrixObservable::~krylovMatrixObservable()
{
	if (dim > 0)
		delete[] tmpBlasVec;
}

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

krylovSpMatrixObservable::krylovSpMatrixObservable(const std::string& name, smatrix* obser) : krylovBasicObservable(name)
{
	dim = obser->m;
	type = SPARSE_MATRIX_TYPE_OBS;
	obs = std::make_unique<smatrix>(*obser); //make owned copy (I don't know if that's the best way)
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
