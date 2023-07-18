# Changelog

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
