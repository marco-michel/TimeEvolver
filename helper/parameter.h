#pragma once

#include <string>
#include <vector>
#include <iomanip>
#include <memory>

/**
* Base class for including model parameters in the output (HDF5) file
*/
class parameter
{
public:
	parameter(const std::string& name, bool inclInFilename) : p_name(name), includeInFilename(inclInFilename) {}
	parameter() = default;
	virtual ~parameter() {}
	virtual std::string getData() = 0;
	virtual bool isDouble() = 0;
	virtual bool isInt() = 0;
	virtual bool isBool() = 0;
	std::string getName() {
		return p_name;
	}
	bool getPrintFilename() {
		return includeInFilename;
	}

private:
	std::string p_name; //!< Name of the model parameter */
	bool includeInFilename; //!< include the parameter also in the filename */
};

/**
* Template class for including model parameters in the output file specialized to parameter datatype
*/
template <typename T>
class typedParameter : public parameter
{
public:
	typedParameter(const std::string& name, bool inclInFilname, const T& data) : parameter(name, inclInFilname), para_data(data) {}

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
#define paraPush(nm, bfnm, ...) std::make_shared<typedParameter<decltype(__VA_ARGS__)>>(nm, bfnm, __VA_ARGS__) //macro for easier insertion