#pragma once

#include <vector>
#include <memory>
#include <string>

#include "parameter.h"
#include "krylovObservables.h"




void saveResult(const std::vector<std::unique_ptr<krylovBasicObservable>>& obs_list, parameter_list& para, const std::string& name);