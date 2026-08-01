#pragma once

#include <vector>
#include <memory>
#include <string>

#include "parameter.h"
#include "krylovObservables.h"
#include "krylovSampleWriter.h"
#include "krylovTimeEvolver.h"


/**
* Auxiliary functions for file output. Some require HDF5.
*
* No HDF5 header is included here; everything HDF5 specific lives in the
* implementation.
*/

/**
* Write the sampled expectation values of a completed time evolution to file.
* The Hamiltonian carried by the result is not written; it can be many
* gigabytes and has its own function, saveSparseMatrix.
* @param result The finished time evolution, providing both the observables and the accuracy information
* @param para List of model parameters, stored as attributes and optionally used for the filename
* @param name Requested filename without extension
*/
void saveResult(const krylovReturn& result, parameter_list& para, const std::string& name);

#ifdef USE_HDF

/**
* The HDF5 result file of one run.
*
* This is what saveResult writes; using the class directly is only necessary in
* order to also record the wavefunction. The sampled states are never held in
* memory, so the file has to be opened before the evolution starts:
*
*     hdf5ResultWriter writer(parameters, "Result");
*     timeEvolver.setSampleWriter(&writer);
*     krylovReturn* result = timeEvolver.timeEvolve();
*     writer.writeObservables(result->observableList);
*     writer.writeMetadata(*result);
*
* Do not call saveResult in addition; it opens the same file for truncation and
* would discard the states that were just written.
*
* The model parameters and the library version are recorded on construction, so
* that an aborted run still leaves behind a file that says what it was.
*/
class hdf5ResultWriter : public TE::krylovSampleWriter
{
public:
    /**
    * Create the result file, replacing an existing one of the same name.
    * @param para List of model parameters, stored as attributes and optionally used for the filename
    * @param name Requested filename without extension
    */
    hdf5ResultWriter(parameter_list& para, const std::string& name);
    ~hdf5ResultWriter() override;

    hdf5ResultWriter(const hdf5ResultWriter&) = delete;
    hdf5ResultWriter& operator=(const hdf5ResultWriter&) = delete;

    //Destination of the sampled wavefunction, called by the time evolution.
    void beginStates(size_t dim, size_t expectedSamples) override;
    void appendState(const std::complex<double>* state, size_t dim) override;
    void finishStates() override;

    /**
    * Write the sampled expectation values, one dataset per observable.
    * @param observables The observables of a finished time evolution
    */
    void writeObservables(const std::vector<std::unique_ptr<krylovBasicObservable>>& observables);

    /**
    * Write what the time evolution reports about itself: the time axis of the
    * samples, the a posteriori error bound, the status code and the size of the
    * computation. Has to be called for the file to be usable, since without the
    * sampling step the samples cannot be placed in time.
    * @param result The finished time evolution
    */
    void writeMetadata(const krylovReturn& result);

    /**
    * The name of the file being written, including the extension and whatever
    * the parameters contributed.
    */
    const std::string& fileName() const;

private:
    struct implementation;
    std::unique_ptr<implementation> impl;
};

/**
* Write a sparse matrix to its own HDF5 file.
* @param mat Matrix to be saved
* @param name Name of the output file
*/
void saveSparseMatrix(const smatrix* mat, const std::string& name);

#endif
