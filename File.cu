#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuComplex.h>

#include "Header.cuh"

// Kernel to compute negativeH = -1.0 * d_HRet_values[j-1+j*m]
// The result is stored in d_negativeH[0]
__global__ void computeNegative(const cuDoubleComplex* d_HRet_values, cuDoubleComplex* d_negativeH, int idx)
{
    // Use a single thread to perform the computation.
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        d_negativeH[0] = cuCmul(make_cuDoubleComplex(-1.0, 0.0), d_HRet_values[idx]);
    }
}



// Exposed interface function using C linkage
extern "C" void launchComputeNegative(const cuDoubleComplex* d_HRet_values,
    cuDoubleComplex* d_negativeH, int idx)
{
    // Launch the kernel
    computeNegative << <1, 1 >> > (d_HRet_values, d_negativeH, idx);
    // Optionally, synchronize to wait for kernel completion
    cudaDeviceSynchronize();
}