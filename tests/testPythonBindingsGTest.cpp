#include <cstdlib>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#ifdef USE_HDF
#include "krylovObservables.h"
#include "krylovTimeEvolver.h"
#include "matrixDataTypes.h"
#endif

#ifndef TE_PYTHON_MODULE_DIR
#error "TE_PYTHON_MODULE_DIR is not defined."
#endif

namespace {

std::string quotedPath(const std::filesystem::path& path) {
    std::string s = path.string();
    if (s.find(' ') == std::string::npos) {
        return s;
    }
    return "'" + s + "'";
}

int runPythonScript(const std::string& scriptBody) {
    const auto scriptPath = std::filesystem::temp_directory_path() / "timeevolver_python_gtest.py";

    std::ofstream scriptFile(scriptPath);
    if (!scriptFile.is_open()) {
        return 1;
    }
    scriptFile << scriptBody;
    scriptFile.close();

    const std::string command = "python3 " + quotedPath(scriptPath);
    const int status = std::system(command.c_str());

    std::error_code ec;
    std::filesystem::remove(scriptPath, ec);
    return status;
}

std::string pythonFloatList(const std::vector<double>& values) {
    std::ostringstream out;
    out << "[";
    out << std::setprecision(17);
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            out << ", ";
        }
        out << values[i];
    }
    out << "]";
    return out.str();
}

}  

TEST(PythonBindings, ModuleImportsAndSymbolsExist) {
    const std::string moduleDir = TE_PYTHON_MODULE_DIR;

    std::string script;
    script += "import sys\n";
    script += "sys.path.insert(0, r'" + moduleDir + "')\n";
    script += "import timeevolver\n";
    script += "assert hasattr(timeevolver, 'TimeEvolver')\n";
    script += "assert hasattr(timeevolver, 'SparseMatrix')\n";
    script += "assert hasattr(timeevolver, 'SpMatrixObservable')\n";
    script += "assert hasattr(timeevolver, 'TimeEvolverResult')\n";
    script += "assert hasattr(timeevolver, 'time_evolve')\n";

    EXPECT_EQ(runPythonScript(script), 0);
}

TEST(PythonBindings, RealTimeEvolutionAndOptionalHDFRoundTrip) {
    const std::string moduleDir = TE_PYTHON_MODULE_DIR;
    const auto tmpDir = std::filesystem::temp_directory_path();
    const auto hdfPath = tmpDir / "timeevolver_python_gtest_matrix.h5";

    std::string script;
    script += "import math\n";
    script += "import os\n";
    script += "import sys\n";
    script += "sys.path.insert(0, r'" + moduleDir + "')\n";
    script += "import timeevolver\n";
    script += "H = timeevolver.SparseMatrix(\n";
    script += "    2, 2,\n";
    script += "    [1+0j, 1+0j],\n";
    script += "    [1, 0],\n";
    script += "    [0, 1],\n";
    script += "    sym=True, hermitian=True\n";
    script += ")\n";
    script += "if hasattr(H, 'save_hdf5'):\n";
    script += "    hdf_file = r'" + hdfPath.string() + "'\n";
    script += "    H.save_hdf5(hdf_file)\n";
    script += "    H = timeevolver.SparseMatrix.load_hdf5(hdf_file)\n";
    script += "    if os.path.exists(hdf_file):\n";
    script += "        os.remove(hdf_file)\n";
    script += "P0 = timeevolver.SparseMatrix(2, 2, [1+0j], [0], [0], hermitian=True)\n";
    script += "obs0 = timeevolver.SpMatrixObservable('P0', P0)\n";
    script += "res = timeevolver.time_evolve(\n";
    script += "    t=math.pi,\n";
    script += "    initial_state=[1+0j, 0+0j],\n";
    script += "    sampling_step=math.pi/2,\n";
    script += "    hamiltonian=H,\n";
    script += "    observables=[obs0],\n";
    script += "    tol=1e-10,\n";
    script += "    m=2,\n";
    script += "    progress_bar=False\n";
    script += ")\n";
    script += "vals = res.observable_values[0]\n";
    script += "assert len(vals) == 3\n";
    script += "expected = [1.0, 0.0, 1.0]\n";
    script += "for v, e in zip(vals, expected):\n";
    script += "    assert abs(v - e) < 1e-6, (v, e)\n";
    script += "assert abs(res.evolved_time - math.pi) < 1e-10\n";
    script += "assert res.status_code in (0, 1, 2)\n";

    EXPECT_EQ(runPythonScript(script), 0);

    std::error_code ec;
    std::filesystem::remove(hdfPath, ec);
}

TEST(PythonBindings, HDF5LoadedMatricesMatchCppCore) {
#ifndef USE_HDF
    GTEST_SKIP() << "HDF5 is not enabled in this build.";
#else
    constexpr double kPi = 3.141592653589793238462643383279502884;
    const std::string moduleDir = TE_PYTHON_MODULE_DIR;
    const auto tmpDir = std::filesystem::temp_directory_path();
    const auto hdfHamiltonianPath = tmpDir / "timeevolver_python_gtest_hamiltonian.h5";
    const auto hdfObservablePath = tmpDir / "timeevolver_python_gtest_observable.h5";

    std::complex<double> hamValues[] = {
        std::complex<double>(1.0, 0.0),
        std::complex<double>(1.0, 0.0),
    };
    size_t hamColumns[] = {1, 0};
    size_t hamRowIndex[] = {0, 1};
    TE::smatrix hamNative(hamValues, hamColumns, hamRowIndex, 2, 2, 2);
    hamNative.sym = true;
    hamNative.hermitian = true;
    hamNative.saveHDF5(hdfHamiltonianPath.string());

    std::complex<double> obsValues[] = {
        std::complex<double>(1.0, 0.0),
    };
    size_t obsColumns[] = {0};
    size_t obsRowIndex[] = {0};
    TE::smatrix obsNative(obsValues, obsColumns, obsRowIndex, 1, 2, 2);
    obsNative.hermitian = true;
    obsNative.saveHDF5(hdfObservablePath.string());

    std::complex<double> initialState[] = {
        std::complex<double>(1.0, 0.0),
        std::complex<double>(0.0, 0.0),
    };

    std::vector<std::unique_ptr<krylovBasicObservable>> nativeObservables;
    nativeObservables.push_back(std::make_unique<krylovSpMatrixObservable>(
        "P0",
        std::make_unique<TE::smatrix>(obsNative)));

    auto hamForEvolution = std::make_unique<TE::smatrix>(hamNative);
    krylovTimeEvolver nativeEvolver(
        kPi,
        initialState,
        kPi / 2.0,
        std::move(nativeObservables),
        std::move(hamForEvolution),
        1.0,
        1.0e-10,
        2,
        false,
        false);

    std::unique_ptr<krylovReturn> nativeResult(nativeEvolver.timeEvolve());
    ASSERT_NE(nativeResult, nullptr);
    ASSERT_EQ(nativeResult->observableList.size(), 1u);

    auto* nativeObservable = nativeResult->observableList[0].get();
    const size_t nativeNumSamples = nativeObservable->retNumSamples();
    std::vector<double> expectedValues(
        nativeObservable->retExpectationValues(),
        nativeObservable->retExpectationValues() + nativeNumSamples);

    std::string script;
    script += "import math\n";
    script += "import sys\n";
    script += "sys.path.insert(0, r'" + moduleDir + "')\n";
    script += "import timeevolver\n";
    script += "H = timeevolver.SparseMatrix.load_hdf5(r'" + hdfHamiltonianPath.string() + "')\n";
    script += "P0 = timeevolver.SparseMatrix.load_hdf5(r'" + hdfObservablePath.string() + "')\n";
    script += "obs0 = timeevolver.SpMatrixObservable('P0', P0)\n";
    script += "res = timeevolver.time_evolve(\n";
    script += "    t=math.pi,\n";
    script += "    initial_state=[1+0j, 0+0j],\n";
    script += "    sampling_step=math.pi/2,\n";
    script += "    hamiltonian=H,\n";
    script += "    observables=[obs0],\n";
    script += "    tol=1e-10,\n";
    script += "    m=2,\n";
    script += "    progress_bar=False\n";
    script += ")\n";
    script += "vals = list(res.observable_values[0])\n";
    script += "expected = " + pythonFloatList(expectedValues) + "\n";
    script += "assert len(vals) == len(expected), (len(vals), len(expected))\n";
    script += "for i, (v, e) in enumerate(zip(vals, expected)):\n";
    script += "    assert abs(v - e) < 1e-8, (i, v, e)\n";

    const int scriptStatus = runPythonScript(script);

    std::error_code ec;
    std::filesystem::remove(hdfHamiltonianPath, ec);
    std::filesystem::remove(hdfObservablePath, ec);

    EXPECT_EQ(scriptStatus, 0);
#endif
}
