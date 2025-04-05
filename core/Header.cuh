#pragma once

#include <cuComplex.h>

#ifdef __cplusplus
extern "C" {
#endif

    // Interface function that launches the kernel.
    void launchComputeNegative(const cuDoubleComplex* d_HRet_values,
        cuDoubleComplex* d_negativeH, int idx);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

    // Interface function that wraps the kernel launch.
    void launchComputeFourOutputsFromDouble(const double* d_input,
        cuDoubleComplex* d_out0,
        cuDoubleComplex* d_out1,
        cuDoubleComplex* d_out2,
        cuDoubleComplex* d_out3);

#ifdef __cplusplus
}
#endif