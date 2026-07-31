#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "krylovObservables.h"
#include "krylovTimeEvolver.h"
#include "matrixDataTypes.h"

namespace py = pybind11;

namespace {

struct PySparseMatrix {
    size_t n = 0;
    size_t m = 0;
    std::vector<std::complex<double>> values;
    std::vector<size_t> columns;
    std::vector<size_t> row_index;
    bool sym = false;
    bool hermitian = false;
    bool upper_tri = false;

    PySparseMatrix() = default;

    PySparseMatrix(
        size_t n_,
        size_t m_,
        std::vector<std::complex<double>> values_,
        std::vector<size_t> columns_,
        std::vector<size_t> row_index_,
        bool sym_ = false,
        bool hermitian_ = false,
        bool upper_tri_ = false)
        : n(n_),
          m(m_),
          values(std::move(values_)),
          columns(std::move(columns_)),
          row_index(std::move(row_index_)),
          sym(sym_),
          hermitian(hermitian_),
          upper_tri(upper_tri_) {
        validate();
    }

    void validate() const {
        if (n == 0 || m == 0) {
            throw py::value_error("SparseMatrix dimensions must be > 0.");
        }
        if (values.empty()) {
            throw py::value_error("SparseMatrix must contain at least one non-zero entry.");
        }
        if (values.size() != columns.size() || values.size() != row_index.size()) {
            throw py::value_error("values, columns, and row_index must have identical lengths.");
        }

        for (size_t i = 0; i < values.size(); ++i) {
            if (columns[i] >= n) {
                throw py::value_error("SparseMatrix column index out of bounds.");
            }
            if (row_index[i] >= m) {
                throw py::value_error("SparseMatrix row_index out of bounds.");
            }
        }
    }

    std::unique_ptr<TE::smatrix> to_native() {
        validate();
        auto mat = std::make_unique<TE::smatrix>(
            values.data(),
            columns.data(),
            row_index.data(),
            values.size(),
            static_cast<unsigned int>(n),
            static_cast<unsigned int>(m));

        mat->sym = sym;
        mat->hermitian = hermitian;
        mat->upperTri = upper_tri;
        return mat;
    }

    void save_hdf5(const std::string& filename) const {
#ifdef USE_HDF
        PySparseMatrix copy = *this;
        auto native = copy.to_native();
        native->saveHDF5(filename);
#else
        (void)filename;
        // pybind11 has no py::runtime_error; std::runtime_error is translated
        // into Python's RuntimeError.
        throw std::runtime_error(
            "HDF5 support is not available in this build (USE_HDF not defined).");
#endif
    }

    static PySparseMatrix load_hdf5(const std::string& filename) {
#ifdef USE_HDF
        TE::smatrix native;
        native.loadHDF5(filename);

        PySparseMatrix out;
        out.n = native.n;
        out.m = native.m;
        out.values.assign(native.values, native.values + native.numValues);
        out.columns.assign(native.columns, native.columns + native.numValues);
        out.row_index.assign(native.rowIndex, native.rowIndex + native.numValues);
        out.sym = native.sym;
        out.hermitian = native.hermitian;
        out.upper_tri = native.upperTri;
        return out;
#else
        (void)filename;
        // pybind11 has no py::runtime_error; std::runtime_error is translated
        // into Python's RuntimeError.
        throw std::runtime_error(
            "HDF5 support is not available in this build (USE_HDF not defined).");
#endif
    }
};

struct PySpMatrixObservable {
    std::string name;
    PySparseMatrix matrix;

    PySpMatrixObservable() = default;

    PySpMatrixObservable(std::string name_, PySparseMatrix matrix_)
        : name(std::move(name_)), matrix(std::move(matrix_)) {
        if (name.empty()) {
            throw py::value_error("Observable name must not be empty.");
        }
    }

    std::unique_ptr<krylovBasicObservable> to_native() {
        return std::make_unique<krylovSpMatrixObservable>(name, matrix.to_native());
    }
};

struct PyResult {
    std::vector<std::complex<double>> evolved_state;
    double err = 0.0;
    double evolved_time = 0.0;
    size_t num_samples = 0;
    size_t n_steps = 0;
    size_t dim = 0;
    size_t krylov_dim = 0;
    int status_code = 0;
    std::vector<std::string> observable_names;
    std::vector<std::vector<double>> observable_values;
};

PyResult convert_result(const krylovReturn& native_result) {
    PyResult result;
    result.err = native_result.err;
    result.evolved_time = native_result.evolvedTime;
    result.num_samples = native_result.numSamples;
    result.n_steps = native_result.n_steps;
    result.dim = native_result.dim;
    result.krylov_dim = native_result.krylovDim;
    result.status_code = native_result.statusCode;

    result.evolved_state.assign(
        native_result.evolvedState,
        native_result.evolvedState + native_result.dim);

    result.observable_names.reserve(native_result.observableList.size());
    result.observable_values.reserve(native_result.observableList.size());

    for (const auto& obs : native_result.observableList) {
        result.observable_names.push_back(obs->retName());
        const double* raw_values = obs->retExpectationValues();
        const size_t n = std::min(obs->retNumSamples(), native_result.numSamples);
        result.observable_values.emplace_back(raw_values, raw_values + n);
    }

    return result;
}

class PyTimeEvolver {
public:
    PyTimeEvolver(
        double t,
        std::vector<std::complex<double>> initial_state,
        double sampling_step,
        PySparseMatrix hamiltonian,
        std::vector<PySpMatrixObservable> observables,
        double exp_factor = 1.0,
        double tol = 1.0e-6,
        int m = 40,
        bool fast_integration = false,
        bool progress_bar = false)
        : t_(t),
          initial_state_(std::move(initial_state)),
          sampling_step_(sampling_step),
          hamiltonian_(std::move(hamiltonian)),
          observables_(std::move(observables)),
          exp_factor_(exp_factor),
          tol_(tol),
          m_(m),
          fast_integration_(fast_integration),
          progress_bar_(progress_bar) {
        validate_inputs();
    }

    PyResult time_evolve() {
        auto native_hamiltonian = hamiltonian_.to_native();

        std::vector<std::unique_ptr<krylovBasicObservable>> native_observables;
        native_observables.reserve(observables_.size());
        for (auto& obs : observables_) {
            native_observables.push_back(obs.to_native());
        }

        krylovTimeEvolver evolver(
            t_,
            initial_state_.data(),
            sampling_step_,
            std::move(native_observables),
            std::move(native_hamiltonian),
            exp_factor_,
            tol_,
            m_,
            fast_integration_,
            progress_bar_);

        std::unique_ptr<krylovReturn> native_result(evolver.timeEvolve());
        return convert_result(*native_result);
    }

private:
    void validate_inputs() const {
        if (t_ <= 0.0) {
            throw py::value_error("t must be > 0.");
        }
        if (sampling_step_ <= 0.0) {
            throw py::value_error("sampling_step must be > 0.");
        }
        if (tol_ <= 0.0) {
            throw py::value_error("tol must be > 0.");
        }
        if (m_ <= 0) {
            throw py::value_error("m must be > 0.");
        }
        if (initial_state_.empty()) {
            throw py::value_error("initial_state must not be empty.");
        }

        hamiltonian_.validate();
        if (hamiltonian_.n != hamiltonian_.m) {
            throw py::value_error("hamiltonian must be square.");
        }
        if (hamiltonian_.m != initial_state_.size()) {
            throw py::value_error(
                "hamiltonian dimension must match len(initial_state).");
        }

        for (const auto& obs : observables_) {
            if (obs.matrix.n != hamiltonian_.n || obs.matrix.m != hamiltonian_.m) {
                throw py::value_error(
                    "Observable dimensions must match hamiltonian dimensions.");
            }
        }

        double norm2 = 0.0;
        for (const auto& value : initial_state_) {
            norm2 += std::norm(value);
        }
        const double norm = std::sqrt(norm2);
        if (std::abs(norm - 1.0) > tol_) {
            throw py::value_error(
                "initial_state must be normalized (within tolerance `tol`).");
        }
    }

    double t_;
    std::vector<std::complex<double>> initial_state_;
    double sampling_step_;
    PySparseMatrix hamiltonian_;
    std::vector<PySpMatrixObservable> observables_;
    double exp_factor_;
    double tol_;
    int m_;
    bool fast_integration_;
    bool progress_bar_;
};

PyResult time_evolve(
    double t,
    std::vector<std::complex<double>> initial_state,
    double sampling_step,
    PySparseMatrix hamiltonian,
    std::vector<PySpMatrixObservable> observables,
    double exp_factor = 1.0,
    double tol = 1.0e-6,
    int m = 40,
    bool fast_integration = false,
    bool progress_bar = false) {
    PyTimeEvolver evolver(
        t,
        std::move(initial_state),
        sampling_step,
        std::move(hamiltonian),
        std::move(observables),
        exp_factor,
        tol,
        m,
        fast_integration,
        progress_bar);
    return evolver.time_evolve();
}

}  // namespace

PYBIND11_MODULE(timeevolver, m) {
    m.doc() = "Python bindings for TimeEvolver evolution";

    py::class_<PySparseMatrix>(m, "SparseMatrix")
        .def(py::init<>())
        .def(py::init<
                 size_t,
                 size_t,
                 std::vector<std::complex<double>>,
                 std::vector<size_t>,
                 std::vector<size_t>,
                 bool,
                 bool,
                 bool>(),
             py::arg("n"),
             py::arg("m"),
             py::arg("values"),
             py::arg("columns"),
             py::arg("row_index"),
             py::arg("sym") = false,
             py::arg("hermitian") = false,
             py::arg("upper_tri") = false)
        .def_readwrite("n", &PySparseMatrix::n)
        .def_readwrite("m", &PySparseMatrix::m)
        .def_readwrite("values", &PySparseMatrix::values)
        .def_readwrite("columns", &PySparseMatrix::columns)
        .def_readwrite("row_index", &PySparseMatrix::row_index)
        .def_readwrite("sym", &PySparseMatrix::sym)
        .def_readwrite("hermitian", &PySparseMatrix::hermitian)
        .def_readwrite("upper_tri", &PySparseMatrix::upper_tri)
        .def("validate", &PySparseMatrix::validate)
        .def("save_hdf5", &PySparseMatrix::save_hdf5, py::arg("filename"))
        .def_static("load_hdf5", &PySparseMatrix::load_hdf5, py::arg("filename"));

    py::class_<PySpMatrixObservable>(m, "SpMatrixObservable")
        .def(py::init<>())
        .def(py::init<std::string, PySparseMatrix>(),
             py::arg("name"),
             py::arg("matrix"))
        .def_readwrite("name", &PySpMatrixObservable::name)
        .def_readwrite("matrix", &PySpMatrixObservable::matrix);

    py::class_<PyResult>(m, "TimeEvolverResult")
        .def_readonly("evolved_state", &PyResult::evolved_state)
        .def_readonly("err", &PyResult::err)
        .def_readonly("evolved_time", &PyResult::evolved_time)
        .def_readonly("num_samples", &PyResult::num_samples)
        .def_readonly("n_steps", &PyResult::n_steps)
        .def_readonly("dim", &PyResult::dim)
        .def_readonly("krylov_dim", &PyResult::krylov_dim)
        .def_readonly("status_code", &PyResult::status_code)
        .def_readonly("observable_names", &PyResult::observable_names)
        .def_readonly("observable_values", &PyResult::observable_values);

    py::class_<PyTimeEvolver>(m, "TimeEvolver")
        .def(py::init<
                 double,
                 std::vector<std::complex<double>>,
                 double,
                 PySparseMatrix,
                 std::vector<PySpMatrixObservable>,
                 double,
                 double,
                 int,
                 bool,
                 bool>(),
             py::arg("t"),
             py::arg("initial_state"),
             py::arg("sampling_step"),
             py::arg("hamiltonian"),
             py::arg("observables"),
             py::arg("exp_factor") = 1.0,
             py::arg("tol") = 1.0e-6,
             py::arg("m") = 40,
             py::arg("fast_integration") = false,
             py::arg("progress_bar") = false)
        .def("time_evolve", &PyTimeEvolver::time_evolve);

    m.def(
        "time_evolve",
        &time_evolve,
        py::arg("t"),
        py::arg("initial_state"),
        py::arg("sampling_step"),
        py::arg("hamiltonian"),
        py::arg("observables"),
        py::arg("exp_factor") = 1.0,
        py::arg("tol") = 1.0e-6,
        py::arg("m") = 40,
        py::arg("fast_integration") = false,
        py::arg("progress_bar") = false,
        "Function to run time evolution.");
}
