#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>

inline const char* cublasStatusToString(cublasStatus_t status)
{
    switch (status) {
    case CUBLAS_STATUS_SUCCESS:
        return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
        return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
        return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
        return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
        return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
        return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
        return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
        return "CUBLAS_STATUS_INTERNAL_ERROR";
    case CUBLAS_STATUS_NOT_SUPPORTED:
        return "CUBLAS_STATUS_NOT_SUPPORTED";
    case CUBLAS_STATUS_LICENSE_ERROR:
        return "CUBLAS_STATUS_LICENSE_ERROR";
    default:
        return "CUBLAS_STATUS_UNKNOWN";
    }
}

inline void throwOnCudaError(cudaError_t status, const char* expr, int line)
{
    if (status == cudaSuccess) {
        return;
    }

    std::ostringstream msg;
    msg << "CUDA API failed at line " << line << " in `" << expr << "` with error: "
        << cudaGetErrorString(status) << " (" << static_cast<int>(status) << ")";
    throw std::runtime_error(msg.str());
}

inline void throwOnCusparseError(cusparseStatus_t status, const char* expr, int line)
{
    if (status == CUSPARSE_STATUS_SUCCESS) {
        return;
    }

    std::ostringstream msg;
    msg << "cuSPARSE API failed at line " << line << " in `" << expr << "` with error: "
        << cusparseGetErrorString(status) << " (" << static_cast<int>(status) << ")";
    throw std::runtime_error(msg.str());
}

inline void throwOnCublasError(cublasStatus_t status, const char* expr, int line)
{
    if (status == CUBLAS_STATUS_SUCCESS) {
        return;
    }

    std::ostringstream msg;
    msg << "cuBLAS API failed at line " << line << " in `" << expr << "` with error: "
        << cublasStatusToString(status) << " (" << static_cast<int>(status) << ")";
    throw std::runtime_error(msg.str());
}

#define TE_CUDA_CHECK(expr) throwOnCudaError((expr), #expr, __LINE__)
#define TE_CUSPARSE_CHECK(expr) throwOnCusparseError((expr), #expr, __LINE__)
#define TE_CUBLAS_CHECK(expr) throwOnCublasError((expr), #expr, __LINE__)
