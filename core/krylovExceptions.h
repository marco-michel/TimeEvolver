#pragma once

#include <stdexcept>
#include <string>

namespace TE {

    /**
    * Base class of every error reported by the TimeEvolver library. Deriving from
    * std::runtime_error means a caller that only wants to know that something went
    * wrong can catch std::exception and read what().
    */
    class krylovError : public std::runtime_error
    {
    public:
        explicit krylovError(const std::string& message) : std::runtime_error(message) {}
    };

    /**
    * A precondition of the called function was violated, e.g. an empty matrix, a
    * state vector that is not normalized or mismatching dimensions. These errors
    * are caused by the input and can be avoided by the caller.
    */
    class krylovInvalidArgument : public krylovError
    {
    public:
        explicit krylovInvalidArgument(const std::string& message) : krylovError(message) {}
    };

    /**
    * The computation could not reach the requested accuracy. This is a legitimate
    * numerical outcome rather than a mistake, so it is worth catching separately:
    * a caller may want to retry with a larger tolerance.
    */
    class krylovConvergenceError : public krylovError
    {
    public:
        explicit krylovConvergenceError(const std::string& message) : krylovError(message) {}
    };

    /**
    * A numerical backend (BLAS, LAPACK or sparse BLAS) reported a failure.
    */
    class krylovBackendError : public krylovError
    {
    public:
        explicit krylovBackendError(const std::string& message) : krylovError(message) {}
    };

}
