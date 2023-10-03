#pragma once

#include <vector>
#include <memory>
#include <string>

#ifdef USE_HDF
#include <H5Cpp.h>
using namespace H5;
#endif

#include "parameter.h"
#include "krylovObservables.h"




/**
* Several auxilary function for file output. Some require HDF5. 
*/

void saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>>& obs_list, parameter_list& para, const std::string& name);

#ifdef USE_HDF
void saveMatrix(const matrix* mat, const std::string& name);
void saveSparseMatrix(const smatrix* mat, const std::string& name);
#endif