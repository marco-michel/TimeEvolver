# Changelog

## 2.0.0 - --

### Added
- Implemented a rudimentary logger (based on BOOST log) including a verbose option. Logging to file will be supported in a later release. 
- To ensure data ownership and avoid the copying of large matrices all large matrices including observables are now expected to be passed as unique_ptr.
- KrylovReturn now includes all output Parameter including the (maybe altered in terms of sorting) Hamiltonian Matrix.
- Intel MKL is now optional (but still highly recommended) and can be replaced by any BLAS library.
- Added wrapper for sparse matrix operations to ensure easier implementation for different sparse BLAS libraries. 
- New class `krylovBasicObservable` for computing and storing expectation values as well as simplyifing output reworking entirely the input/output of the TimeEvolver.
- Progressbar is now in a seperate thread in order not to reduce computing performance.
- `krylovBasicObservable` can throw `requestStopException` exception to stop the time evolution.
- Experimental Mac support.
- Version header `version.h`.
- Model parameters `typedParameter` have an additional bool variable to indicate if its value should also be included in the output filename.
  
### Changed
- Imput parameter `expFactor` of TimeEvolver is now a real (not complex) datatype.
- Reworked constructor of `TimeEvolver` with API-breaking effects. 
- Expectation value computation is now outsourced from the TimeEvolver core to `krylovBasicObservable`.
- C++ version requirement changed to 17.
- Initial vector is required to be normalized to 1. 

### Fixed
- Fixed some memory leaks
- Added checks for empty matrices
- Fixed a confusing output behavior where `dim` in `krylovReturn` stores not the Hilbertspace dim suggested by its name but the Krylov dimension.
- Implemented exception handling for failed numerical integration 

## 1.4.1 - 2022-04-11

### Fixed
- Fixed memory leak and memory corruption in smatrix class

### Changed
- Code layout for observable class

## 1.4.0 - 2022-03-07

### Added
- New `krylovBasicObservable` class to extend types of observables (vector, sparseMatrix, denseMatrix)
- `krylovTimeEvolver` can relay computation of expectation value to observable `krylovBasicObservable` instance 
- Occupation number in `basis` can be restricted, e.g. only even numbers
- Progressbar

### Changed
- No memory reservation in basis creation


## 1.3.0 - 2022-10-25

### Added
- New `outputHelper` class to simplify output in h5 or csv files

### Changed
- Rework of `createMatrix` for better performance and less memory consumption
- Minimal C++ std increased to C++14 for compatibility with future boost release

## 1.2.0 - 2022-10-11

### Added
- Added ctest to check compilation

## 1.1.1 - 2022-09-19

### Fixed
- Changed access specifier of the method `createMatrix` from private to protected


## 1.1.0 - 2022-03-30

### Added
- Added a new example: Two coupled oscillators
- Added installation routine in CMake files

### Fixed
- Made HDF5 optional
