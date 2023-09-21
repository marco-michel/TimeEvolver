#pragma once

#include <sstream>
#include <iomanip>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>


/**
* Logger used to print warnings, errors, etc.
*/
class krylovLogger
{
public:
	enum loggingLevel { DEBUG, INFO, WARNING, ERROR, FATAL };

	krylovLogger() = default;
	krylovLogger(loggingLevel level);
	void set_loggingLevel(loggingLevel level);
	void log_message(loggingLevel level, const std::string& msg);
	static std::string stringStreamWrapperd(double val);

private:
	bool logToFile;
	std::string logFileName;
	boost::log::trivial::severity_level translateLevel(loggingLevel level);

};