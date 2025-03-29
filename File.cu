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

__global__ void computeFourOutputsFromDouble(const double* input,
    cuDoubleComplex* out0,
    cuDoubleComplex* out1,
    cuDoubleComplex* out2,
    cuDoubleComplex* out3)
{
    if (threadIdx.x == 0 && blockIdx.x == 0)
    {
        double d = *input;  // Read the input value.

        // Convert input to a complex number.
        cuDoubleComplex cVal = make_cuDoubleComplex(d, 0.0);

        // Write the two outputs that are just the input value.
        *out0 = cVal;
        *out1 = cVal;

        // Compute and write the inverse: 1/d.
        *out2 = make_cuDoubleComplex(1.0 / d, 0.0);

        // Compute and write the negative: -d.
        *out3 = make_cuDoubleComplex(-d, 0.0);
    }
}

// Wrapper function callable from standard C++ code.
extern "C" void launchComputeFourOutputsFromDouble(const double* d_input,
    cuDoubleComplex* d_out0,
    cuDoubleComplex* d_out1,
    cuDoubleComplex* d_out2,
    cuDoubleComplex* d_out3)
{
    // Launch the kernel with one block and one thread.
    computeFourOutputsFromDouble << <1, 1 >> > (d_input, d_out0, d_out1, d_out2, d_out3);
    cudaDeviceSynchronize();
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