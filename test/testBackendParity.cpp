#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "Basis.h"
#include "exampleHamiltonian.h"
#include "hamiltonian.h"
#include "krylovObservables.h"
#include "krylovTimeEvolver.h"
#include "matrixDataTypes.h"

using namespace TE;

namespace {

struct parityCase {
    std::unique_ptr<smatrix> hamiltonian;
    std::vector<std::unique_ptr<smatrix>> observableMatrices;
    std::vector<std::complex<double>> initialState;
    double maxT = 0.0;
    double samplingStep = 0.0;
    double tol = 0.0;
    int krylovDim = 0;
    bool fastIntegration = false;
    double expFactor = 1.0;
};

std::string cudaSkipReason()
{
#ifdef USE_CUDA
    int deviceCount = 0;
    const cudaError_t countStatus = cudaGetDeviceCount(&deviceCount);
    if (countStatus != cudaSuccess) {
        return std::string("CUDA unavailable: ") + cudaGetErrorString(countStatus);
    }
    if (deviceCount == 0) {
        return "CUDA unavailable: no device detected";
    }
    const cudaError_t setStatus = cudaSetDevice(0);
    if (setStatus != cudaSuccess) {
        return std::string("CUDA unavailable: ") + cudaGetErrorString(setStatus);
    }
    return std::string();
#else
    return "CUDA support is not compiled in";
#endif
}

std::vector<std::complex<double>> createDeterministicVector(std::size_t length)
{
    std::vector<std::complex<double>> values(length);
    for (std::size_t i = 0; i < length; ++i) {
        const double phase = static_cast<double>(i + 1);
        values[i] = std::complex<double>(std::sin(0.13 * phase), std::cos(0.19 * phase));
    }

    const double norm = cblas_dznrm2(static_cast<int>(length), values.data(), 1);
    const std::complex<double> invNorm(1.0 / norm, 0.0);
    cblas_zscal(static_cast<int>(length), &invNorm, values.data(), 1);
    return values;
}

double maxVectorDifference(const std::complex<double>* lhs, const std::complex<double>* rhs, std::size_t length)
{
    double maxDiff = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        maxDiff = std::max(maxDiff, std::abs(lhs[i] - rhs[i]));
    }
    return maxDiff;
}

double maxVectorDifference(const std::vector<std::complex<double>>& lhs, const std::vector<std::complex<double>>& rhs)
{
    if (lhs.size() != rhs.size()) {
        return std::numeric_limits<double>::infinity();
    }
    return maxVectorDifference(lhs.data(), rhs.data(), lhs.size());
}

double maxScalarDifference(const double* lhs, const double* rhs, std::size_t length)
{
    double maxDiff = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        maxDiff = std::max(maxDiff, std::abs(lhs[i] - rhs[i]));
    }
    return maxDiff;
}

std::vector<std::unique_ptr<krylovBasicObservable>> cloneObservables(const std::vector<std::unique_ptr<smatrix>>& observableMatrices)
{
    std::vector<std::unique_ptr<krylovBasicObservable>> observables;
    observables.reserve(observableMatrices.size());
    for (std::size_t i = 0; i < observableMatrices.size(); ++i) {
        observables.push_back(std::make_unique<krylovSpMatrixObservable>(
            "mode" + std::to_string(i),
            std::make_unique<smatrix>(*observableMatrices[i])));
    }
    return observables;
}

parityCase buildSimpleParityCase()
{
    const int particleCount = 40;
    const int numberModes = 2;
    const double gap1 = 1.0;
    const double gap2 = 1.7;
    const double coupling = 0.4;

    parityCase result;
    result.maxT = 2.0;
    result.samplingStep = 0.25;
    result.tol = 1.0e-10;
    result.krylovDim = 20;
    result.fastIntegration = false;

    basis simpleBasis(particleCount, numberModes, 0, 0);
    Hamiltonian hamiltonian;
    hamiltonian.hamiltonOperator.push_back(hamiltonian.createNumberOperator(0, gap1));
    hamiltonian.hamiltonOperator.push_back(hamiltonian.createNumberOperator(1, gap2));
    hamiltonian.hamiltonOperator.push_back(hamiltonian.linInteraction(0, 1, 0, 0, false, coupling));
    hamiltonian.hamiltonOperator.push_back(hamiltonian.linInteraction(0, 1, 0, 0, true, coupling));

    result.hamiltonian = hamiltonian.createHamiltonMatrix(&simpleBasis);
    result.observableMatrices = hamiltonian.createNumberOperatorObservables(&simpleBasis);
    result.initialState.assign(simpleBasis.numberElements, std::complex<double>(0.0, 0.0));

    basisVector init(numberModes);
    init.e[0] = particleCount;
    const int entry = simpleBasis.hashTable.find(init)->second;
    result.initialState[entry] = std::complex<double>(1.0, 0.0);

    return result;
}

parityCase buildBlackHoleParityCase()
{
    const int N0 = 8;
    const int Nm = 2;
    const int K = 4;

    parityCase result;
    result.maxT = 1.5;
    result.samplingStep = 0.3;
    result.tol = 1.0e-9;
    result.krylovDim = 16;
    result.fastIntegration = false;

    tensorBasis basisStateSpace(N0, 2, Nm, 2 * K, 1);
    exampleHamiltonian hamiltonian(N0, Nm, 4.0, K, K, 1, 0.03, 1.0, 1.0, 0.08, 0.08);
    hamiltonian.createSimplifiedHamiltonian();

    result.hamiltonian = hamiltonian.createHamiltonMatrix(&basisStateSpace);
    result.observableMatrices = hamiltonian.createNumberOperatorObservables(&basisStateSpace);
    result.initialState.assign(basisStateSpace.numberElements, std::complex<double>(0.0, 0.0));

    basisVector init = hamiltonian.createInitState();
    const int entry = basisStateSpace.hashTable.find(init)->second;
    result.initialState[entry] = std::complex<double>(1.0, 0.0);

    return result;
}

std::unique_ptr<krylovReturn> runTimeEvolution(const parityCase& input, krylovTimeEvolver::executionBackend backend)
{
    std::vector<std::complex<double>> initialState = input.initialState;
    auto observables = cloneObservables(input.observableMatrices);

    krylovTimeEvolver evolver(
        input.maxT,
        initialState.data(),
        input.samplingStep,
        std::move(observables),
        std::make_unique<smatrix>(*input.hamiltonian),
        input.expFactor,
        input.tol,
        input.krylovDim,
        input.fastIntegration,
        false,
        backend);
    evolver.changeLogLevel(krylovLogger::FATAL);
    return std::unique_ptr<krylovReturn>(evolver.timeEvolve());
}

void compareReturns(const krylovReturn& cpu, const krylovReturn& gpu, double stateTolerance, double observableTolerance)
{
    ASSERT_EQ(cpu.statusCode, gpu.statusCode);
    ASSERT_EQ(cpu.numSamples, gpu.numSamples);
    ASSERT_EQ(cpu.observableList.size(), gpu.observableList.size());
    EXPECT_NEAR(cpu.evolvedTime, gpu.evolvedTime, 1.0e-12);
    EXPECT_NEAR(cpu.err, gpu.err, std::max(1.0e-12, 1.0e-6 * std::abs(cpu.err)));

    const double stateDiff = maxVectorDifference(cpu.evolvedState, gpu.evolvedState, cpu.dim);
    EXPECT_LE(stateDiff, stateTolerance) << "Evolved state mismatch.";

    for (std::size_t i = 0; i < cpu.observableList.size(); ++i) {
        SCOPED_TRACE(cpu.observableList[i]->retName());
        ASSERT_EQ(cpu.observableList[i]->retNumSamples(), gpu.observableList[i]->retNumSamples());
        const double observableDiff = maxScalarDifference(
            cpu.observableList[i]->retExpectationValues(),
            gpu.observableList[i]->retExpectationValues(),
            cpu.observableList[i]->retNumSamples());
        EXPECT_LE(observableDiff, observableTolerance) << "Observable mismatch.";
    }
}

class BackendParityTest : public ::testing::Test {};

TEST_F(BackendParityTest, SpMVMatchesCPU)
{
    const std::string reason = cudaSkipReason();
    if (!reason.empty()) {
        GTEST_SKIP() << reason;
    }

    parityCase input = buildBlackHoleParityCase();
    std::vector<std::complex<double>> x = createDeterministicVector(input.hamiltonian->m);
    std::vector<std::complex<double>> cpuOut(input.hamiltonian->n, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> gpuOut(input.hamiltonian->n, std::complex<double>(0.0, 0.0));

    input.hamiltonian->spMV(smatrix::one, x.data(), cpuOut.data());

#ifdef USE_CUDA
    smatrixCUDA gpuMatrix(*input.hamiltonian);
    const std::size_t xBytes = x.size() * sizeof(cuDoubleComplex);
    const std::size_t yBytes = gpuOut.size() * sizeof(cuDoubleComplex);
    TE_CUDA_CHECK(cudaMemcpy(gpuMatrix.CX, x.data(), xBytes, cudaMemcpyHostToDevice));
    gpuMatrix.spMV(gpuMatrix.oneCUDA, gpuMatrix.vecX, gpuMatrix.vecY);
    TE_CUDA_CHECK(cudaDeviceSynchronize());
    TE_CUDA_CHECK(cudaMemcpy(gpuOut.data(), gpuMatrix.CY, yBytes, cudaMemcpyDeviceToHost));
#endif

    const double diff = maxVectorDifference(cpuOut, gpuOut);
    EXPECT_LE(diff, 1.0e-11);
}

TEST_F(BackendParityTest, SimpleTimeEvolutionMatchesCPU)
{
    const std::string reason = cudaSkipReason();
    if (!reason.empty()) {
        GTEST_SKIP() << reason;
    }

    const parityCase input = buildSimpleParityCase();
    const std::unique_ptr<krylovReturn> cpu = runTimeEvolution(input, krylovTimeEvolver::executionBackend::CPU);
    const std::unique_ptr<krylovReturn> gpu = runTimeEvolution(input, krylovTimeEvolver::executionBackend::GPU);
    compareReturns(*cpu, *gpu, 1.0e-10, 1.0e-10);
}

TEST_F(BackendParityTest, BlackHoleTimeEvolutionMatchesCPU)
{
    const std::string reason = cudaSkipReason();
    if (!reason.empty()) {
        GTEST_SKIP() << reason;
    }

    const parityCase input = buildBlackHoleParityCase();
    const std::unique_ptr<krylovReturn> cpu = runTimeEvolution(input, krylovTimeEvolver::executionBackend::CPU);
    const std::unique_ptr<krylovReturn> gpu = runTimeEvolution(input, krylovTimeEvolver::executionBackend::GPU);
    compareReturns(*cpu, *gpu, 1.0e-8, 1.0e-8);
}

} // namespace
