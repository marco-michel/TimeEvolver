#pragma once

#include <vector>
#include <memory>
#include <iomanip>

#include "krylovTimeEvolver.h"
//#include "krylovHelper.h"

#ifdef USE_HDF
#include <H5Cpp.h>
using namespace H5;
#endif



//------------------------------------------------------------------------------------------------

class observable  //intended to be extened to something similar as the parameter list
{
public:
	observable(const std::string& o_name, obsType o_type): name(o_name), type(o_type) {}
	std::string getName()
	{
		return name;
	}
	obsType getType()
	{
		return type;
	}
private:
	obsType type;
	std::string name;
};

typedef std::vector< std::shared_ptr<observable> > observable_list;
#define obsPush(nm, ...) std::make_shared<observable>(nm, __VA_ARGS__) //macro for easier insertion

//------------------------------------------------------------------------------------------------


//Helper to outsource writing results to file
class outputHelper
{
public: 
	outputHelper(krylovReturn* res, const parameter_list& para, const observable_list& obs, const std::string& name) : results(res), para_list(para), obs_list(obs), fileName(name) { nbObservables = obs_list.size(); nbOutputParameters = para_list.size(); }
	void saveResult();

private:
	krylovReturn* results; parameter_list para_list; observable_list obs_list; std::string fileName;  
	unsigned int nbObservables; unsigned int nbOutputParameters;
#ifdef USE_HDF
	DataSet dataset;
	void writeAttributes();
#endif
};


