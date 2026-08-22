#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "Basis.h"
#include "exampleHamiltonian.h"
#include "hamiltonian.h"
#include "krylovLogger.h"
#include "krylovObservables.h"
#include "krylovTimeEvolver.h"
#include "matrixDataTypes.h"

#ifdef USE_MKL
#include <mkl.h>
#endif

namespace po = boost::program_options;
using namespace TE;

namespace {

enum class benchmarkKind { spmv, arnoldi, e2e };
enum class backendMode { cpu, gpu, both };

struct blackHoleCaseParameters {
    std::string name;
    int N0 = 20;
    int Nm = 4;
    int K = 8;
    double DeltaN = 12.0;
    int capacity = 1;
    double C0 = 0.0003;
    double Cm = 0.065;
    double maxT = 1000.0;
    double samplingStep = 20.0;
    double tol = 1.0e-6;
    int krylovDim = 40;
    bool fastIntegration = true;
};

struct benchmarkCase {
    std::string name;
    bool hermitianStorage = false;
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

struct timingSummary {
    double totalMs = 0.0;
    double minMs = 0.0;
    double maxMs = 0.0;
    double meanMs = 0.0;
    double medianMs = 0.0;
};

struct benchmarkRecord {
    std::string benchmark;
    std::string caseName;
    std::string backend;
    std::size_t dimension = 0;
    std::size_t nnz = 0;
    double avgNnzPerRow = std::numeric_limits<double>::quiet_NaN();
    double matrixHostGiB = std::numeric_limits<double>::quiet_NaN();
    double matrixDeviceGiB = std::numeric_limits<double>::quiet_NaN();
    double krylovBasisGiB = std::numeric_limits<double>::quiet_NaN();
    int krylovDim = 0;
    int warmup = 0;
    int repeats = 0;
    int observables = 0;
    timingSummary wall;
    double deviceMeanMs = std::numeric_limits<double>::quiet_NaN();
    double transferMeanMs = std::numeric_limits<double>::quiet_NaN();
    double maxError = std::numeric_limits<double>::quiet_NaN();
    std::string note;
};

void populateStaticMetrics(benchmarkRecord& record);

benchmarkRecord makeDescriptionRecord(const benchmarkCase& input, const std::string& backendLabel)
{
    benchmarkRecord description;
    description.benchmark = "case";
    description.caseName = input.name;
    description.backend = backendLabel;
    description.dimension = input.hamiltonian->n;
    description.nnz = input.hamiltonian->numValues;
    description.krylovDim = input.krylovDim;
    description.observables = static_cast<int>(input.observableMatrices.size());
    populateStaticMetrics(description);
    return description;
}

std::string benchmarkToString(benchmarkKind benchmark)
{
    switch (benchmark) {
    case benchmarkKind::spmv:
        return "spmv";
    case benchmarkKind::arnoldi:
        return "arnoldi";
    case benchmarkKind::e2e:
        return "e2e";
    }

    throw std::invalid_argument("Unknown benchmark kind.");
}

benchmarkKind parseBenchmarkKind(const std::string& value)
{
    if (value == "spmv") {
        return benchmarkKind::spmv;
    }
    if (value == "arnoldi") {
        return benchmarkKind::arnoldi;
    }
    if (value == "e2e") {
        return benchmarkKind::e2e;
    }

    throw std::invalid_argument("Unsupported benchmark `" + value + "`.");
}

backendMode parseBackendMode(const std::string& value)
{
    if (value == "cpu") {
        return backendMode::cpu;
    }
    if (value == "gpu") {
        return backendMode::gpu;
    }
    if (value == "both") {
        return backendMode::both;
    }

    throw std::invalid_argument("Unsupported backend `" + value + "`.");
}

double bytesToGiB(long double bytes)
{
    static constexpr long double gib = 1024.0L * 1024.0L * 1024.0L;
    return static_cast<double>(bytes / gib);
}

double estimateSparseMatrixHostGiB(std::size_t nnz)
{
    const long double bytesPerValue =
        sizeof(std::complex<double>) + sizeof(std::size_t) + sizeof(std::size_t);
    return bytesToGiB(static_cast<long double>(nnz) * bytesPerValue);
}

double estimateSparseMatrixDeviceGiB(std::size_t rows, std::size_t nnz)
{
    const long double bytesPerValue =
        sizeof(std::complex<double>) + sizeof(std::size_t);
    const long double rowOffsetBytes =
        static_cast<long double>(rows + 1U) * sizeof(std::size_t);
    return bytesToGiB(static_cast<long double>(nnz) * bytesPerValue + rowOffsetBytes);
}

double estimateKrylovBasisGiB(std::size_t dimension, int krylovDim)
{
    const long double bytes =
        static_cast<long double>(dimension) *
        static_cast<long double>(krylovDim) *
        sizeof(std::complex<double>);
    return bytesToGiB(bytes);
}

void populateStaticMetrics(benchmarkRecord& record)
{
    if (record.dimension != 0U) {
        record.avgNnzPerRow = static_cast<double>(record.nnz) / static_cast<double>(record.dimension);
    }
    record.matrixHostGiB = estimateSparseMatrixHostGiB(record.nnz);
    record.matrixDeviceGiB = estimateSparseMatrixDeviceGiB(record.dimension, record.nnz);
    record.krylovBasisGiB = estimateKrylovBasisGiB(record.dimension, record.krylovDim);
}

timingSummary summarizeSamples(const std::vector<double>& samples)
{
    if (samples.empty()) {
        throw std::invalid_argument("No timing samples collected.");
    }

    timingSummary summary;
    summary.totalMs = std::accumulate(samples.begin(), samples.end(), 0.0);
    summary.minMs = *std::min_element(samples.begin(), samples.end());
    summary.maxMs = *std::max_element(samples.begin(), samples.end());
    summary.meanMs = summary.totalMs / static_cast<double>(samples.size());

    std::vector<double> sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    const std::size_t mid = sorted.size() / 2;
    if ((sorted.size() % 2U) == 0U) {
        summary.medianMs = 0.5 * (sorted[mid - 1] + sorted[mid]);
    } else {
        summary.medianMs = sorted[mid];
    }

    return summary;
}

std::vector<std::complex<double>> createDenseBenchmarkVector(std::size_t length)
{
    std::vector<std::complex<double>> values(length);

    for (std::size_t i = 0; i < length; ++i) {
        const double phase = static_cast<double>(i + 1);
        values[i] = std::complex<double>(std::sin(0.17 * phase), std::cos(0.11 * phase));
    }

    const double norm = cblas_dznrm2(static_cast<int>(length), values.data(), 1);
    const std::complex<double> invNorm(1.0 / norm, 0.0);
    cblas_zscal(static_cast<int>(length), &invNorm, values.data(), 1);
    return values;
}

benchmarkCase buildSimpleCase(bool hermitianStorage)
{
    const int particleCount = 200;
    const int numberModes = 2;
    const double gap1 = 1.0;
    const double gap2 = 2.0;
    const double coupling = 1.0;

    benchmarkCase result;
    result.name = "simple";
    result.maxT = 10.0;
    result.samplingStep = 0.01;
    result.tol = 1.0e-6;
    result.krylovDim = 40;
    result.fastIntegration = false;
    result.expFactor = 1.0;

    basis simpleBasis(particleCount, numberModes, 0, 0);
    Hamiltonian hamiltonian;
    std::vector<opTerm> hamiltonianTerms;
    hamiltonianTerms.push_back(hamiltonian.createNumberOperator(0, gap1));
    hamiltonianTerms.push_back(hamiltonian.createNumberOperator(1, gap2));
    hamiltonianTerms.push_back(hamiltonian.linInteraction(0, 1, 0, 0, false, coupling));
    hamiltonianTerms.push_back(hamiltonian.linInteraction(0, 1, 0, 0, true, coupling));
    hamiltonian.hamiltonOperator = hamiltonianTerms;

    result.hamiltonian = hamiltonian.createHamiltonMatrix(&simpleBasis, hermitianStorage);
    result.observableMatrices = hamiltonian.createNumberOperatorObservables(&simpleBasis);
    result.initialState.assign(simpleBasis.numberElements, std::complex<double>(0.0, 0.0));

    basisVector init(numberModes);
    init.e[0] = particleCount;
    const int entry = simpleBasis.hashTable.find(init)->second;
    result.initialState[entry] = std::complex<double>(1.0, 0.0);

    return result;
}

std::string defaultBlackHoleCaseName(const blackHoleCaseParameters& parameters)
{
    return "blackhole_N0" + std::to_string(parameters.N0) +
        "_Nm" + std::to_string(parameters.Nm) +
        "_K" + std::to_string(parameters.K);
}

blackHoleCaseParameters makePresetBlackHoleCase(
    const std::string& name,
    int N0,
    int Nm,
    int K,
    double C0,
    double Cm,
    double maxT,
    double samplingStep,
    double tol,
    int krylovDim,
    bool fastIntegration)
{
    blackHoleCaseParameters parameters;
    parameters.name = name;
    parameters.N0 = N0;
    parameters.Nm = Nm;
    parameters.K = K;
    parameters.C0 = C0;
    parameters.Cm = Cm;
    parameters.maxT = maxT;
    parameters.samplingStep = samplingStep;
    parameters.tol = tol;
    parameters.krylovDim = krylovDim;
    parameters.fastIntegration = fastIntegration;
    return parameters;
}

benchmarkCase buildBlackHoleCase(const blackHoleCaseParameters& parameters, bool buildObservables, bool hermitianStorage)
{
    benchmarkCase result;
    result.name = parameters.name.empty() ? defaultBlackHoleCaseName(parameters) : parameters.name;
    result.maxT = parameters.maxT;
    result.samplingStep = parameters.samplingStep;
    result.tol = parameters.tol;
    result.krylovDim = parameters.krylovDim;
    result.fastIntegration = parameters.fastIntegration;
    result.expFactor = 1.0;

    tensorBasis basisStateSpace(parameters.N0, 2, parameters.Nm, 2 * parameters.K, parameters.capacity);
    exampleHamiltonian hamiltonian(
        parameters.N0,
        parameters.Nm,
        parameters.DeltaN,
        parameters.K,
        parameters.K,
        parameters.capacity,
        parameters.C0,
        1.0,
        1.0,
        parameters.Cm,
        parameters.Cm);
    hamiltonian.createSimplifiedHamiltonian();

    result.hamiltonian = hamiltonian.createHamiltonMatrix(&basisStateSpace, hermitianStorage);
    if (buildObservables) {
        result.observableMatrices = hamiltonian.createNumberOperatorObservables(&basisStateSpace);
    }
    result.initialState.assign(basisStateSpace.numberElements, std::complex<double>(0.0, 0.0));

    basisVector init = hamiltonian.createInitState();
    const int entry = basisStateSpace.hashTable.find(init)->second;
    result.initialState[entry] = std::complex<double>(1.0, 0.0);

    return result;
}

benchmarkCase buildCaseByName(const std::string& caseName, bool buildObservables, const blackHoleCaseParameters& customParameters, bool hermitianStorage)
{
    if (caseName == "simple") {
        return buildSimpleCase(hermitianStorage);
    }
    if (caseName == "blackhole_lb") {
        return buildBlackHoleCase(makePresetBlackHoleCase(caseName, 1, 1, 1, 1.0, 1.0, 1000.0, 1.0, 1.0e-8, 40, false), buildObservables, hermitianStorage);
    }
    if (caseName == "blackhole_medium") {
        return buildBlackHoleCase(makePresetBlackHoleCase(caseName, 20, 2, 4, 1.0, 1.0, 10.0, 0.01, 1.0e-8, 40, false), buildObservables, hermitianStorage);
    }
    if (caseName == "blackhole_large") {
        return buildBlackHoleCase(makePresetBlackHoleCase(caseName, 20, 4, 8, 0.0003, 0.065, 1000.0, 20.0, 1.0e-6, 40, true), buildObservables, hermitianStorage);
    }
    if (caseName == "blackhole_quick") {
        return buildBlackHoleCase(makePresetBlackHoleCase(caseName, 20, 5, 9, 0.0003, 0.065, 10.0, 1.0, 1.0e-6, 40, true), buildObservables, hermitianStorage);
    }
    if (caseName == "blackhole_1g") {
        return buildBlackHoleCase(makePresetBlackHoleCase(caseName, 24, 6, 9, 0.0003, 0.065, 10.0, 1.0, 1.0e-6, 40, true), buildObservables, hermitianStorage);
    }
    if (caseName == "blackhole_2g") {
        return buildBlackHoleCase(makePresetBlackHoleCase(caseName, 20, 6, 10, 0.0003, 0.065, 40.0, 1.0, 1.0e-6, 40, true), buildObservables, hermitianStorage);
    }
    if (caseName == "custom") {
        return buildBlackHoleCase(customParameters, buildObservables, hermitianStorage);
    }

    throw std::invalid_argument("Unknown benchmark case `" + caseName + "`.");
}

std::vector<std::string> availableCaseNames()
{
    return { "simple", "blackhole_lb", "blackhole_medium", "blackhole_large", "blackhole_quick", "blackhole_1g", "blackhole_2g", "custom" };
}

bool optionWasProvided(const po::variables_map& vm, const char* optionName)
{
    const auto found = vm.find(optionName);
    return found != vm.end() && !found->second.defaulted();
}

void applyRuntimeOverrides(
    benchmarkCase& input,
    const blackHoleCaseParameters& overrides,
    const po::variables_map& vm)
{
    if (optionWasProvided(vm, "maxT")) {
        input.maxT = overrides.maxT;
    }
    if (optionWasProvided(vm, "samplingStep")) {
        input.samplingStep = overrides.samplingStep;
    }
    if (optionWasProvided(vm, "tol")) {
        input.tol = overrides.tol;
    }
    if (optionWasProvided(vm, "m")) {
        input.krylovDim = overrides.krylovDim;
    }
    if (optionWasProvided(vm, "fastIntegration")) {
        input.fastIntegration = overrides.fastIntegration;
    }
}

std::vector<std::unique_ptr<krylovBasicObservable>> cloneObservables(const benchmarkCase& input, bool includeObservables)
{
    std::vector<std::unique_ptr<krylovBasicObservable>> observables;
    if (!includeObservables) {
        return observables;
    }

    observables.reserve(input.observableMatrices.size());
    for (std::size_t i = 0; i < input.observableMatrices.size(); ++i) {
        observables.push_back(std::make_unique<krylovSpMatrixObservable>(
            "mode" + std::to_string(i),
            std::make_unique<smatrix>(*input.observableMatrices[i])));
    }

    return observables;
}

template <typename Callable>
timingSummary measureWallClock(int warmup, int repeats, Callable&& callable)
{
    for (int i = 0; i < warmup; ++i) {
        callable();
    }

    std::vector<double> samples;
    samples.reserve(repeats);

    for (int i = 0; i < repeats; ++i) {
        const auto start = std::chrono::steady_clock::now();
        callable();
        const auto stop = std::chrono::steady_clock::now();
        const std::chrono::duration<double, std::milli> elapsed = stop - start;
        samples.push_back(elapsed.count());
    }

    return summarizeSamples(samples);
}


benchmarkRecord runCpuSpmvBenchmark(const benchmarkCase& input, int warmup, int repeats)
{
    std::vector<std::complex<double>> hostInput = createDenseBenchmarkVector(input.hamiltonian->m);
    std::vector<std::complex<double>> hostOutput(input.hamiltonian->n, std::complex<double>(0.0, 0.0));

    benchmarkRecord result;
    result.benchmark = "spmv";
    result.caseName = input.name;
    result.backend = "cpu";
    result.dimension = input.hamiltonian->n;
    result.nnz = input.hamiltonian->numValues;
    result.krylovDim = input.krylovDim;
    result.warmup = warmup;
    result.repeats = repeats;
    populateStaticMetrics(result);

    result.wall = measureWallClock(warmup, repeats, [&]() {
        input.hamiltonian->spMV(smatrix::one, hostInput.data(), hostOutput.data());
    });

    return result;
}

class arnoldiBenchmarkHarness : public krylovTimeEvolver
{
public:
    explicit arnoldiBenchmarkHarness(const benchmarkCase& input)
        : krylovTimeEvolver(
            input.maxT,
            const_cast<std::complex<double>*>(input.initialState.data()),
            input.samplingStep,
            {},
            std::make_unique<smatrix>(*input.hamiltonian),
            input.expFactor,
            input.tol,
            input.krylovDim,
            input.fastIntegration,
            false)
    {
        changeLogLevel(krylovLogger::FATAL);
    }

    void runArnoldiStep()
    {
        TE::matrix hessenberg(m, m);
        TE::matrix basis(Hsize, m);
        double h = 0.0;
        std::size_t actualKrylovDim = 0;
        const double tolRate = tol / t;

        arnoldiAlgorithm(tolRate, &hessenberg, &basis, &h, &actualKrylovDim);
    }
};

benchmarkRecord runArnoldiBenchmark(const benchmarkCase& input, int warmup, int repeats, const std::string& backendLabel)
{
    arnoldiBenchmarkHarness harness(input);

    benchmarkRecord result;
    result.benchmark = "arnoldi";
    result.caseName = input.name;
    result.backend = backendLabel;
    result.dimension = input.hamiltonian->n;
    result.nnz = input.hamiltonian->numValues;
    result.krylovDim = input.krylovDim;
    result.warmup = warmup;
    result.repeats = repeats;
    populateStaticMetrics(result);

    result.wall = measureWallClock(warmup, repeats, [&]() {
        harness.runArnoldiStep();
    });

    return result;
}

benchmarkRecord runE2EBenchmark(const benchmarkCase& input, int warmup, int repeats, const std::string& backendLabel, bool includeObservables)
{
    if (input.samplingStep > input.maxT) {
        throw std::invalid_argument(
            "Invalid e2e benchmark case `" + input.name +
            "`: samplingStep > maxT, so timeEvolve() would only sample the initial state and skip the evolution loop.");
    }

    benchmarkRecord result;
    result.benchmark = "e2e";
    result.caseName = input.name;
    result.backend = backendLabel;
    result.dimension = input.hamiltonian->n;
    result.nnz = input.hamiltonian->numValues;
    result.krylovDim = input.krylovDim;
    result.warmup = warmup;
    result.repeats = repeats;
    result.observables = includeObservables ? static_cast<int>(input.observableMatrices.size()) : 0;
    populateStaticMetrics(result);

    result.wall = measureWallClock(warmup, repeats, [&]() {
        std::vector<std::complex<double>> initialState = input.initialState;
        std::unique_ptr<smatrix> hamiltonian = std::make_unique<smatrix>(*input.hamiltonian);
        std::vector<std::unique_ptr<krylovBasicObservable>> observables = cloneObservables(input, includeObservables);

        krylovTimeEvolver evolver(
            input.maxT,
            initialState.data(),
            input.samplingStep,
            std::move(observables),
            std::move(hamiltonian),
            input.expFactor,
            input.tol,
            input.krylovDim,
            input.fastIntegration,
            false);
        evolver.changeLogLevel(krylovLogger::FATAL);

        std::unique_ptr<krylovReturn> output(evolver.timeEvolve());
    });

    return result;
}

void printRecord(const benchmarkRecord& record)
{
    std::cout << std::fixed << std::setprecision(3);
    std::cout << record.benchmark << " case=" << record.caseName
              << " backend=" << record.backend
              << " dim=" << record.dimension
              << " nnz=" << record.nnz
              << " nnz_per_row=" << record.avgNnzPerRow
              << " matrix_host_gib=" << record.matrixHostGiB
              << " matrix_device_gib=" << record.matrixDeviceGiB
              << " basis_gib=" << record.krylovBasisGiB
              << " m=" << record.krylovDim
              << " repeats=" << record.repeats
              << " warmup=" << record.warmup
              << " mean_ms=" << record.wall.meanMs
              << " median_ms=" << record.wall.medianMs
              << " min_ms=" << record.wall.minMs
              << " max_ms=" << record.wall.maxMs;

    if (!std::isnan(record.deviceMeanMs)) {
        std::cout << " device_mean_ms=" << record.deviceMeanMs;
    }
    if (!std::isnan(record.transferMeanMs)) {
        std::cout << " transfer_mean_ms=" << record.transferMeanMs;
    }
    if (!std::isnan(record.maxError)) {
        std::cout << " max_error=" << record.maxError;
    }
    if (!record.note.empty()) {
        std::cout << " note=\"" << record.note << "\"";
    }
    std::cout << std::endl;
}

void appendCsvRecord(const benchmarkRecord& record, const std::string& csvPath)
{
    if (csvPath.empty()) {
        return;
    }

    bool writeHeader = false;
    {
        std::ifstream existing(csvPath);
        writeHeader = !existing.good() || existing.peek() == std::ifstream::traits_type::eof();
    }

    std::ofstream output(csvPath, std::ios::app);
    if (!output.is_open()) {
        throw std::runtime_error("Failed to open CSV output `" + csvPath + "`.");
    }

    if (writeHeader) {
        output << "benchmark,case,backend,dimension,nnz,avg_nnz_per_row,matrix_host_gib,matrix_device_gib,krylov_basis_gib,krylov_dim,warmup,repeats,observables,total_ms,mean_ms,median_ms,min_ms,max_ms,device_mean_ms,transfer_mean_ms,max_error,note\n";
    }

    output << record.benchmark << ','
           << record.caseName << ','
           << record.backend << ','
           << record.dimension << ','
           << record.nnz << ','
           << record.avgNnzPerRow << ','
           << record.matrixHostGiB << ','
           << record.matrixDeviceGiB << ','
           << record.krylovBasisGiB << ','
           << record.krylovDim << ','
           << record.warmup << ','
           << record.repeats << ','
           << record.observables << ','
           << record.wall.totalMs << ','
           << record.wall.meanMs << ','
           << record.wall.medianMs << ','
           << record.wall.minMs << ','
           << record.wall.maxMs << ','
           << record.deviceMeanMs << ','
           << record.transferMeanMs << ','
           << record.maxError << ','
           << '"' << record.note << '"' << '\n';
}

std::vector<benchmarkRecord> describeCase(const benchmarkCase& input, backendMode backend)
{
    std::vector<benchmarkRecord> records;

    if (backend == backendMode::cpu || backend == backendMode::both) {
        benchmarkRecord cpuRecord = makeDescriptionRecord(input, "cpu");
        cpuRecord.note = "matrix assembly is host-side";
        records.push_back(std::move(cpuRecord));
    }

    if (backend == backendMode::gpu || backend == backendMode::both) {
        if (backend == backendMode::gpu) {
            throw std::runtime_error("GPU case description requested but CUDA support is not compiled in.");
        }

        benchmarkRecord skipped = makeDescriptionRecord(input, "gpu");
        skipped.note = "skipped: CUDA support is not compiled in";
        records.push_back(std::move(skipped));
    }

    return records;
}

std::vector<benchmarkRecord> runBenchmarks(const benchmarkCase& input, benchmarkKind benchmark, backendMode backend, int warmup, int repeats, bool includeObservables)
{
    std::vector<benchmarkRecord> records;

    const auto maybeRunCpu = [&]() {
        switch (benchmark) {
        case benchmarkKind::spmv:
            records.push_back(runCpuSpmvBenchmark(input, warmup, repeats));
            return;
        case benchmarkKind::arnoldi:
            records.push_back(runArnoldiBenchmark(input, warmup, repeats, "cpu"));
            return;
        case benchmarkKind::e2e:
            records.push_back(runE2EBenchmark(input, warmup, repeats, "cpu", includeObservables));
            return;
        }
    };

    const auto maybeRunGpu = [&]() {
        throw std::runtime_error("GPU benchmarks requested but CUDA support is not compiled in.");
    };

    if (backend == backendMode::cpu || backend == backendMode::both) {
        maybeRunCpu();
    }
    if (backend == backendMode::gpu) {
        maybeRunGpu();
    }
    if (backend == backendMode::both) {
        try {
            maybeRunGpu();
        }
        catch (const std::exception& e)
        {
            benchmarkRecord skipped;
            skipped.benchmark = benchmarkToString(benchmark);
            skipped.caseName = input.name;
            skipped.backend = "gpu";
            skipped.dimension = input.hamiltonian->n;
            skipped.nnz = input.hamiltonian->numValues;
            skipped.krylovDim = input.krylovDim;
            skipped.warmup = warmup;
            skipped.repeats = repeats;
            skipped.observables = includeObservables ? static_cast<int>(input.observableMatrices.size()) : 0;
            populateStaticMetrics(skipped);
            skipped.note = std::string("skipped: ") + e.what();
            records.push_back(skipped);
        }
    }

    return records;
}

} // namespace

int main(int argc, char* argv[])
{
    std::string benchmarkArg;
    std::string backendArg;
    std::string caseArg;
    std::string csvPath;
    blackHoleCaseParameters customParameters;
    int repeats = 10;
    int warmup = 2;
    int cpuThreads = 0;
    bool hermitianStorage = false;
    bool includeObservables = false;
    bool listCases = false;
    bool describeCaseFlag = false;

    po::options_description options("Benchmark options");
    options.add_options()
        ("help", "Show benchmark options")
        ("benchmark", po::value<std::string>(&benchmarkArg)->default_value("spmv"), "Benchmark to run: spmv, arnoldi, e2e")
        ("backend", po::value<std::string>(&backendArg)->default_value("cpu"), "Backend to run: cpu, gpu, both")
        ("case", po::value<std::string>(&caseArg)->default_value("simple"), "Benchmark case: simple, blackhole_lb, blackhole_medium, blackhole_large, blackhole_quick, blackhole_1g, blackhole_2g, custom")
        ("repeats", po::value<int>(&repeats)->default_value(10), "Measured repetitions")
        ("warmup", po::value<int>(&warmup)->default_value(2), "Warmup repetitions")
        ("csv", po::value<std::string>(&csvPath)->default_value(""), "Append CSV output to this file")
        ("cpu-threads", po::value<int>(&cpuThreads)->default_value(0), "Set MKL CPU thread count if available")
        ("hermitian-storage", po::value<bool>(&hermitianStorage)->default_value(false), "Store only the upper triangle of the Hamiltonian")
        ("include-observables", po::value<bool>(&includeObservables)->default_value(false), "Include observable evaluation in e2e benchmark")
        ("list-cases", po::bool_switch(&listCases), "List available benchmark cases")
        ("describe-case", po::bool_switch(&describeCaseFlag), "Build the case and print dimension, nnz, and memory estimates without timing")
        ("N0", po::value<int>(&customParameters.N0)->default_value(customParameters.N0), "Custom case: number of particles in control sector")
        ("Nm", po::value<int>(&customParameters.Nm)->default_value(customParameters.Nm), "Custom case: number of particles in critical sector")
        ("K", po::value<int>(&customParameters.K)->default_value(customParameters.K), "Custom case: number of critical modes per sector")
        ("DeltaN", po::value<double>(&customParameters.DeltaN)->default_value(customParameters.DeltaN), "Custom case: distance between critical sectors")
        ("capacity", po::value<int>(&customParameters.capacity)->default_value(customParameters.capacity), "Custom case: maximal occupation of critical modes")
        ("C0", po::value<double>(&customParameters.C0)->default_value(customParameters.C0), "Custom case: coupling in control sector")
        ("Cm", po::value<double>(&customParameters.Cm)->default_value(customParameters.Cm), "Custom case: coupling in critical sector")
        ("maxT", po::value<double>(&customParameters.maxT)->default_value(customParameters.maxT), "Simulation time; overrides named blackhole presets too")
        ("samplingStep", po::value<double>(&customParameters.samplingStep)->default_value(customParameters.samplingStep), "Time interval of sampling; overrides named blackhole presets too")
        ("tol", po::value<double>(&customParameters.tol)->default_value(customParameters.tol), "Numerical tolerance; overrides named blackhole presets too")
        ("m", po::value<int>(&customParameters.krylovDim)->default_value(customParameters.krylovDim), "Krylov-space dimension; overrides named blackhole presets too")
        ("fastIntegration", po::value<bool>(&customParameters.fastIntegration)->default_value(customParameters.fastIntegration), "Use faster and less accurate integration; overrides named blackhole presets too");

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, options), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << options << std::endl;
        return 1;
    }

    if (vm.count("help")) {
        std::cout << options << std::endl;
        return 0;
    }

    if (listCases) {
        for (const std::string& availableCase : availableCaseNames()) {
            std::cout << availableCase << std::endl;
        }
        return 0;
    }

    if (repeats <= 0 || warmup < 0) {
        std::cerr << "Warmup must be >= 0 and repeats must be > 0." << std::endl;
        return 1;
    }

#ifdef USE_MKL
    if (cpuThreads > 0) {
        mkl_set_num_threads(cpuThreads);
    }
#else
    (void) cpuThreads;
#endif

    try {
        const benchmarkKind benchmark = parseBenchmarkKind(benchmarkArg);
        const backendMode backend = parseBackendMode(backendArg);
        if (customParameters.name.empty()) {
            customParameters.name = defaultBlackHoleCaseName(customParameters);
        }

        const bool buildObservables = includeObservables && benchmark == benchmarkKind::e2e;
        benchmarkCase selectedCase = buildCaseByName(caseArg, buildObservables, customParameters, hermitianStorage);
        applyRuntimeOverrides(selectedCase, customParameters, vm);

        if (describeCaseFlag) {
            const std::vector<benchmarkRecord> descriptions = describeCase(selectedCase, backend);
            for (const benchmarkRecord& description : descriptions) {
                printRecord(description);
                appendCsvRecord(description, csvPath);
            }
            return 0;
        }

        const std::vector<benchmarkRecord> records = runBenchmarks(selectedCase, benchmark, backend, warmup, repeats, includeObservables);

        for (const benchmarkRecord& record : records) {
            printRecord(record);
            appendCsvRecord(record, csvPath);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Benchmark failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
