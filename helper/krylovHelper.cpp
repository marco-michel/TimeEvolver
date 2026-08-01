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
    * Root attributes the writer itself uses. A model parameter of the same name
    * would be rejected by HDF5 with an error that says nothing about the cause,
    * so the collision is reported here instead. samplingStep and totalTime are
    * absent because those are replaced rather than refused, see writeMetadata.
    */
    const char* const reservedNames[] = { "timeEvolverVersion", "timeEvolverVersionNumber",
        "err", "evolvedTime", "statusCode", "nSteps", "krylovDim", "dim", "numSamples" };

    bool isReserved(const std::string& name)
    {
        for (const char* reserved : reservedNames)
        {
            if (name == reserved)
                return true;
        }
        return false;
    }

    /**
    * Attach one model parameter as a scalar attribute.
    */
    void writeParameter(const H5::H5Object& object, const std::shared_ptr<parameter>& para)
    {
        if (isReserved(para->getName()))
            throw TE::krylovIOError("Parameter '" + para->getName() +
                "' cannot be written, because the time evolver records an attribute of that name itself. Please rename it.");

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

    //The wavefunction is written straight out of the sampled array, without an
    //intermediate buffer, which requires the layouts to agree.
    static_assert(sizeof(TE::hdf5::complexType) == sizeof(std::complex<double>),
        "The complex compound type has to match the memory layout of std::complex<double>");

    const char* const wavefunctionName = "wavefunction";

#endif

}

#ifdef USE_HDF

/**
* Everything that requires an HDF5 type, kept out of the header so that
* including krylovHelper.h does not require HDF5 to be present.
*/
struct hdf5ResultWriter::implementation
{
    std::string fileName;
    H5::H5File file;

    H5::DataSet states;
    H5::CompType complexType;
    //Length of one state. Zero while no wavefunction is being written.
    hsize_t dim = 0;
    //Number of states written so far, which is also the size of the dataset.
    hsize_t rows = 0;
    bool statesOpen = false;

    implementation(const std::string& name)
        : fileName(name), file(name, H5F_ACC_TRUNC), complexType(TE::hdf5::complexDataType()) {}
};

/**
* Create the result file and record what is already known about the run
* @param para List of model parameters
* @param name Requested filename without extension
*/
hdf5ResultWriter::hdf5ResultWriter(parameter_list& para, const std::string& name)
{
    std::string outputFileName = buildFileName(para, name) + ".h5";

TE_HDF5_TRY
    impl.reset(new implementation(outputFileName));

    TE::hdf5::writeVersion(impl->file);
    for (parameter_list::iterator iter = para.begin(); iter != para.end(); iter++)
        writeParameter(impl->file, *iter);
TE_HDF5_CATCH("Could not create " + outputFileName)
}

hdf5ResultWriter::~hdf5ResultWriter() = default;

const std::string& hdf5ResultWriter::fileName() const
{
    return impl->fileName;
}

/**
* Create the dataset that receives the wavefunction. It is extendible rather
* than of fixed size, because an evolution stopped by an observable delivers
* fewer samples than were planned.
* @param dim Number of components of a single state
* @param expectedSamples Number of samples the time evolution intends to take
*/
void hdf5ResultWriter::beginStates(size_t dim, size_t expectedSamples)
{
    if (dim == 0)
        throw TE::krylovInvalidArgument("hdf5ResultWriter: a state cannot have zero components.");

TE_HDF5_TRY
    hsize_t columns = (hsize_t)dim;
    hsize_t initial[2] = { 0, columns };
    hsize_t maximal[2] = { H5S_UNLIMITED, columns };

    H5::DataSpace space(2, initial, maximal);
    impl->states = impl->file.createDataSet(wavefunctionName, impl->complexType, space,
        TE::hdf5::extendibleDatasetProperties((hsize_t)expectedSamples, columns, sizeof(TE::hdf5::complexType)),
        TE::hdf5::extendibleDatasetAccess());

    TE::hdf5::writeAttribute(impl->states, "dim", dim);

    impl->dim = columns;
    impl->rows = 0;
    impl->statesOpen = true;
TE_HDF5_CATCH("Could not create the wavefunction dataset in " + impl->fileName)
}

/**
* Append one sampled state as the next row of the wavefunction dataset
* @param state The state at the current sampling point
* @param dim Number of components
*/
void hdf5ResultWriter::appendState(const std::complex<double>* state, size_t dim)
{
    if (!impl->statesOpen)
        throw TE::krylovIOError("hdf5ResultWriter: a state was offered before the wavefunction dataset was opened.");
    if ((hsize_t)dim != impl->dim)
        throw TE::krylovInvalidArgument("hdf5ResultWriter: the length of the state changed during the time evolution.");

TE_HDF5_TRY
    hsize_t extended[2] = { impl->rows + 1, impl->dim };
    impl->states.extend(extended);

    hsize_t count[2] = { 1, impl->dim };
    hsize_t offset[2] = { impl->rows, 0 };

    H5::DataSpace fileSpace = impl->states.getSpace();
    fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
    H5::DataSpace memorySpace(2, count);

    impl->states.write(state, impl->complexType, memorySpace, fileSpace);
    impl->rows++;
TE_HDF5_CATCH("Could not write sample " + std::to_string(impl->rows) + " to " + impl->fileName)
}

/**
* Close the wavefunction dataset and record how many states it ended up holding
*/
void hdf5ResultWriter::finishStates()
{
    if (!impl->statesOpen)
        return;

TE_HDF5_TRY
    TE::hdf5::writeAttribute(impl->states, "numSamples", (size_t)impl->rows);
    impl->states.close();
    impl->statesOpen = false;

    //Make sure the states are on disk and readable even if the process does not
    //survive to write the metadata.
    impl->file.flush(H5F_SCOPE_GLOBAL);
TE_HDF5_CATCH("Could not close the wavefunction dataset in " + impl->fileName)
}

/**
* Write the sampled expectation values, one dataset per observable
* @param observables The observables of a finished time evolution
*/
void hdf5ResultWriter::writeObservables(const std::vector<std::unique_ptr<krylovBasicObservable>>& observables)
{
TE_HDF5_TRY
    for (auto obsIter = observables.begin(); obsIter != observables.end(); obsIter++)
    {
        hsize_t numSamples = (hsize_t)(*obsIter)->retNumSamples();

        H5::DataSpace space(1, &numSamples);
        H5::DataSet dataset = impl->file.createDataSet((*obsIter)->retName(),
            H5::PredType::IEEE_F64LE, space,
            TE::hdf5::datasetProperties(numSamples, sizeof(double)));

        if (numSamples != 0)
            dataset.write((*obsIter)->retExpectationValues(), H5::PredType::NATIVE_DOUBLE);

        //The run parameters are attributes of the file, not of each observable.
        TE::hdf5::writeAttribute(dataset, "numSamples", (*obsIter)->retNumSamples());
        TE::hdf5::writeAttribute(dataset, "observableType", (int)(*obsIter)->retType());
    }
TE_HDF5_CATCH("Could not write the observables to " + impl->fileName)
}

/**
* Write what the time evolution reports about itself. All of it belongs to the
* file as a whole rather than to any single dataset.
* @param result The finished time evolution
*/
void hdf5ResultWriter::writeMetadata(const krylovReturn& result)
{
TE_HDF5_TRY
    TE::hdf5::writeAttribute(impl->file, "err", result.err);
    TE::hdf5::writeAttribute(impl->file, "evolvedTime", result.evolvedTime);

    //The time axis of the samples, which are evenly spaced. A model parameter
    //of the same name is replaced, since the value the time evolution used is
    //the authoritative one.
    TE::hdf5::writeOrReplaceAttribute(impl->file, "samplingStep", result.samplingStep);
    TE::hdf5::writeOrReplaceAttribute(impl->file, "totalTime", result.totalTime);
    TE::hdf5::writeAttribute(impl->file, "statusCode", result.statusCode);
    TE::hdf5::writeAttribute(impl->file, "nSteps", result.n_steps);
    TE::hdf5::writeAttribute(impl->file, "krylovDim", result.krylovDim);
    TE::hdf5::writeAttribute(impl->file, "dim", result.dim);
    TE::hdf5::writeAttribute(impl->file, "numSamples", result.numSamples);
TE_HDF5_CATCH("Could not write the run metadata to " + impl->fileName)
}

#endif

/**
* Write the sampled expectation values of a completed time evolution to file.
* @param result The finished time evolution
* @param para List of model parameters
* @param name Requested filename without extension
*/
void saveResult(const krylovReturn& result, parameter_list& para, const std::string& name)
{
#ifdef USE_HDF

    hdf5ResultWriter writer(para, name);
    writer.writeObservables(result.observableList);
    writer.writeMetadata(result);

    //If HDF5 is not available, write data to simple csv files
#else

    std::string outputFileName = buildFileName(para, name);
    const std::vector<std::unique_ptr<krylovBasicObservable>>& obs_list = result.observableList;

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
