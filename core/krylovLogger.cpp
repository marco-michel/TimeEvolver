#include "krylovLogger.h"

/**
* Constructor for the logger. 
* @param level Initial logging level. This sets the threshold on the severity for which messages will be printed on the screen. 
*/
krylovLogger::krylovLogger(loggingLevel level)
{
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= translateLevel(level));
}


/**
* Function to change the logging level.
* @param level Logging level. This sets the threshold on the severity for which messages will be printed on the screen.
*/
void krylovLogger::set_loggingLevel(loggingLevel level)
{
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= translateLevel(level));
}

/**
* Function to log messages.
* @param level Severity level of the message, e.g. debuginfo, warning, error, ...
* @param msg Actual message
*/
void krylovLogger::log_message(loggingLevel level, const std::string& msg)
{
    switch (level)
    {
    case DEBUG:     BOOST_LOG_TRIVIAL(debug) << msg; break;
    case INFO:      BOOST_LOG_TRIVIAL(info) << msg;  break;
    case WARNING:   BOOST_LOG_TRIVIAL(warning) << msg;  break;
    case ERROR:     BOOST_LOG_TRIVIAL(error) << msg; break;
    case FATAL:     BOOST_LOG_TRIVIAL(fatal) << msg; 
    }
}


/**
* Internal function to translate logging level to boost logging level
* @param level TimeEvolver logging level
* @return Boost logging level
*/
boost::log::trivial::severity_level krylovLogger::translateLevel(loggingLevel level)
{
    switch (level) {
    case DEBUG:     return boost::log::trivial::severity_level::debug;
    case INFO:      return boost::log::trivial::severity_level::info;
    case WARNING:   return boost::log::trivial::severity_level::warning;
    case ERROR:     return boost::log::trivial::severity_level::error;
    case FATAL:     return boost::log::trivial::severity_level::fatal;
    default:        return boost::log::trivial::severity_level::fatal;
    }
}

/**
* Wrapper for creating single line stringstream outputs
* @param val Number to be converted into a string
* @return string representation of the passed double 
*/
std::string krylovLogger::stringStreamWrapperd(double val)
{
    std::stringstream ss{};
    ss << val;
    return ss.str();
}
