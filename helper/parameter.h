#pragma once

#include <string>
#include <vector>
#include <iomanip>
#include <memory>

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
	std::string getName()
	{
		return p_name;
	}

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