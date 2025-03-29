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