#include <fstream>

#include "krylovHelper.h"

#ifdef USE_HDF
#include "hdf5Support.h"
#endif

namespace {

    /**
    * Build the output filename from the requested name and those parameters that
    * asked to appear in it.
    */
    std::string buildFileName(parameter_list& para, const std::string& name)
    {
        std::string fileName = name;
        for (parameter_list::iterator iter = para.begin(); iter != para.end(); iter++)
        {
            if ((*iter)->getPrintFilename() == true)
                fileName += "_" + (*iter)->getName() + (*iter)->getData();
        }
        return fileName;
    }

#ifdef USE_HDF

    /**
    * Attach one model parameter as a scalar attribute.
    */
    void writeParameter(const H5::H5Object& object, const std::shared_ptr<parameter>& para)
    {
        if (para->isDouble())
            TE::hdf5::writeAttribute(object, para->getName(), dynamic_cast<typedParameter<double>&>(*para).getValue());
        else if (para->isInt())
            TE::hdf5::writeAttribute(object, para->getName(), dynamic_cast<typedParameter<int>&>(*para).getValue());
        else if (para->isBool())
            TE::hdf5::writeAttribute(object, para->getName(), dynamic_cast<typedParameter<bool>&>(*para).getValue());
        else
            throw TE::krylovIOError("Parameter '" + para->getName() +
                "' has a type that cannot be written to file. Supported types are double, int and bool.");
    }

    /**
    * Attach everything describing the run: the accuracy information from the
    * time evolution, the model parameters and the library version. All of it
    * belongs to the file as a whole, not to any single observable.
    */
    void writeRunMetadata(const H5::H5File& file, const krylovReturn& result, parameter_list& para)
    {
        TE::hdf5::writeVersion(file);

        //The a posteriori error bound is the point of the method, so it is part
        //of the result rather than something the caller has to remember to keep.
        TE::hdf5::writeAttribute(file, "err", result.err);
        TE::hdf5::writeAttribute(file, "evolvedTime", result.evolvedTime);
        TE::hdf5::writeAttribute(file, "statusCode", result.statusCode);
        TE::hdf5::writeAttribute(file, "nSteps", result.n_steps);
        TE::hdf5::writeAttribute(file, "krylovDim", result.krylovDim);
        TE::hdf5::writeAttribute(file, "dim", result.dim);
        TE::hdf5::writeAttribute(file, "numSamples", result.numSamples);

        for (parameter_list::iterator iter = para.begin(); iter != para.end(); iter++)
            writeParameter(file, *iter);
    }

#endif

}

/**
* Write the sampled expectation values of a completed time evolution to file.
* @param result The finished time evolution
* @param para List of model parameters
* @param name Requested filename without extension
*/
void saveResult(const krylovReturn& result, parameter_list& para, const std::string& name)
{
    std::string outputFileName = buildFileName(para, name);
    const std::vector<std::unique_ptr<krylovBasicObservable>>& obs_list = result.observableList;

#ifdef USE_HDF

    outputFileName += ".h5";

TE_HDF5_TRY
    H5::H5File file(outputFileName, H5F_ACC_TRUNC);

    writeRunMetadata(file, result, para);

    for (auto obsIter = obs_list.begin(); obsIter != obs_list.end(); obsIter++)
    {
        hsize_t numSamples = (hsize_t)(*obsIter)->retNumSamples();

        H5::DataSpace space(1, &numSamples);
        H5::DataSet dataset = file.createDataSet((*obsIter)->retName(),
            H5::PredType::IEEE_F64LE, space,
            TE::hdf5::datasetProperties(numSamples, sizeof(double)));

        if (numSamples != 0)
            dataset.write((*obsIter)->retExpectationValues(), H5::PredType::NATIVE_DOUBLE);

        //Only what genuinely differs between observables belongs here; the run
        //parameters are attributes of the file.
        TE::hdf5::writeAttribute(dataset, "numSamples", (*obsIter)->retNumSamples());
        TE::hdf5::writeAttribute(dataset, "observableType", (int)(*obsIter)->retType());
    }
TE_HDF5_CATCH("Could not write result to " + outputFileName)

    //If HDF5 is not available, write data to simple csv files
#else

    for (auto obsIter = obs_list.begin(); obsIter != obs_list.end(); obsIter++)
    {
        std::string fileNameCSV = outputFileName + (*obsIter)->retName() + ".csv";
        std::ofstream outputfile;
        outputfile.open(fileNameCSV);

        const double* expectationValues = (*obsIter)->retExpectationValues();
        const size_t numSamples = (*obsIter)->retNumSamples();
        if (numSamples != 0)
        {
            for (size_t i = 0; i + 1 != numSamples; i++)
                outputfile << expectationValues[i] << ", ";
            outputfile << expectationValues[numSamples - 1];
        }
        outputfile.close();
    }

#endif

}


#ifdef USE_HDF
/**
* Write a sparse matrix to its own HDF5 file. The format is the one smatrix
* reads back, so that anything written here can be loaded again.
* @param mat Matrix to be saved
* @param name Name of the output file
*/
void saveSparseMatrix(const smatrix* mat, const std::string& name)
{
    if (mat == nullptr)
        throw TE::krylovInvalidArgument("saveSparseMatrix: no matrix given.");

    mat->saveHDF5(name);
}
#endif
