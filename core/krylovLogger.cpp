#include "krylovLogger.h"


krylovLogger::krylovLogger(loggingLevel level)
{
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= translateLevel(level));
    logToFile = false;
}

void krylovLogger::set_loggingLevel(loggingLevel level)
{
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= translateLevel(level));
}

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
* @param level logging level
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
*/
std::string krylovLogger::stringStreamWrapperd(double val)
{
    std::stringstream ss{};
    ss << val;
    return ss.str();
}
