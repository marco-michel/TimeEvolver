#pragma once

#include <vector>
#include <memory>
#include <string>

#include "parameter.h"
#include "krylovObservables.h"
#include "krylovTimeEvolver.h"


/**
* Auxiliary functions for file output. Some require HDF5.
*
* Note that no HDF5 header is included here: everything HDF5 specific lives in
* the implementation, so that including this header does not drag H5 into the
* including translation unit.
*/

/**
* Write the sampled expectation values of a completed time evolution to file.
* The Hamiltonian carried by the result is deliberately not written; it can be
* many gigabytes and has its own function, saveSparseMatrix.
* @param result The finished time evolution, providing both the observables and the accuracy information
* @param para List of model parameters, stored as attributes and optionally used for the filename
* @param name Requested filename without extension
*/
void saveResult(const krylovReturn& result, parameter_list& para, const std::string& name);

#ifdef USE_HDF
/**
* Write a sparse matrix to its own HDF5 file.
* @param mat Matrix to be saved
* @param name Name of the output file
*/
void saveSparseMatrix(const smatrix* mat, const std::string& name);
#endif
