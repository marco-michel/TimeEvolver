#pragma once

#include <vector>
#include <memory>
#include <iomanip>

#include "krylovTimeEvolver.h"
#include "krylovHelper.h"

#ifdef USE_HDF
#include <H5Cpp.h>
using namespace H5;
#endif


//parent class
class parameter
{
public:
	parameter(const std::string& name) : p_name(name) {}
	virtual ~parameter() {}
	virtual std::string getData() = 0;
	virtual bool isDouble() = 0;
	virtual bool isInt() = 0;
	virtual bool isBool() = 0;
	std::string getName();

private:
	std::string p_name;
};

//parameter class for different data types
template <typename T>
class typedParameter : public parameter
{
public:
	typedParameter(const std::string& name, const T& data) : parameter(name), para_data(data) {}

	std::string getData() //write to string
	{
		std::ostringstream oss;
		oss << std::setprecision(10) << para_data; //set amount of digits in the output string
		return oss.str();
	}

	bool isDouble()
	{
		return std::is_same<T, double>::value;
	}
	
	bool isInt()
	{
		return std::is_same<T, int>::value;
	}
	
	bool isBool()
	{
		return std::is_same<T, bool>::value;
	}

	
	T getValue() 
	{
		return para_data;
	}

private:
	T para_data;
};


typedef std::vector< std::shared_ptr<parameter> > parameter_list; //shorten notation.

#define paraPush(nm, ...) std::make_shared<typedParameter<decltype(__VA_ARGS__)>>(nm, __VA_ARGS__) //macro for easier insertion

//------------------------------------------------------------------------------------------------

enum obsType {VECTOR_TYPE_OBS, SPARSE_MATRIX_TYPE_OBS, MATRIX_TYPE_OBS}; 
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


