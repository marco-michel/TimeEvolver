# Introduction


This is the TimeEvolver software package. 

# Project Structure

The TimeEvolver distribution includes the following files and directories:

<pre>
README                          this file
LICENSE                         the MIT licence
core                            the methods for time evolution using Krylov subspace techniques
example                         concrete example to demonstrate usage of the program
helper                          methods to create the set of basis states and to compute the Hamiltonian matrix
</pre>


# External Requirements

This software packages relies on following external libraries:
* Intel Math Kernel Library (MKL) (required)
* HDF5 (optional)
* BOOST (optional)

On top of that it uses the build managing software Cmake.

## Linux 

We recommend using Linux, or a virtual Linux machine, for TimeEvolver. The above libraries and software packages are usually available via a package manager like ``apt``. For example on Ubuntu can be installed via 
```
sudo apt-get install git
```
```
sudo apt-get install cmake
```
```
sudo apt-get install intel-mkl-full libgsl-dev libhdf5-dev libboost-program-options-dev 
```
Note: If you install Intel MKL via the package manager you will be asked if you want to make MKL your default BLAS/LAPACK library. That is not necessary, so you can choose the default answer "No". Additionally, older versions of ``cmake`` might not find the oneapi version of mkl. We therefore recommend cmake version 3.15 or newer. 

## Mac

To compile the source code of the TimeEvolver the command line developer tools are required. The installation can be prompted by typing e.g. ``g++`` in a terminal. After that you can install the Intel MKL library via the oneAPI installer provided by Intel. The HDF5 library can either be installed with the pre-compiled package (We recommend this option, however, it requires a free registration) or compiled from source (cmake version). In case the library is installed in a custome location please set the variable ``HDF_DIR5`` in the CMakeLists.txt file to the path containing the ``hdf5-targets.cmake`` file. A example of this is already included as a comment in the file.
The boost library can be compiled from source following the instructions found on the corresponding website. To compile and install the library system-wide one can use
```
cd path/to/boost_1_76_0
./bootstrap.sh --prefix=/usr/local
./b2
sudo ./b2 install
```
We note that cmake was not able to locate the mkl library properly via the  command ``find_packge(BLAS)`` within our setting. A workaround is to use the ``FindMKL.cmake`` package that is shipped with the standalone version of mkl. The required folder ``examples_core_c.tgz`` can be found (zipped) in the examples folder in the oneapi directory, which is usually located at ``/opt/intel/oneapi/mkl/latest/example``. Unzip this archive and copy the folder ``examples_core_c/cmake`` into the main TimeEvolver directory. Furthermore, the following adjustments to CMakeLists.txt are required in that case: 
* Add: list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
* Add: include(check_param)
* Replace: find_package(BLAS REQUIRED) -> find_package(MKL MODULE REQUIRED)
* Set: BLAS_LIBRARIES with MKL_LINK_PREFIX

We provide a commented out template in the CMakeLists.txt file. Note also the comments on the standalone version of Intel oneapi below. 

## Windows 

BOOST and HDF5 can be obtained for Windows 64bit OS for example with the ``vcpkg`` package manager via
```
.\vcpkg.exe install hdf5:x64-windows
.\vcpkg.exe install boost:x64-windows
```
However, manual adjustments in cmake may be required. We found it easier to fall back to a virtual linux machine with minuscule performance loss.  

## Standalone Version of Intel MKL

If you use the standalone version of the MKL library, which is now part of the oneAPI framework, please make sure that you follow all installation steps including setting environment variables. For root installation this can be done with:
```
source /opt/intel/oneapi/setvars.sh intel64
```
Note that the variables are only set for the context of your session. For a permanent solution please visit the Intel helppage. 


# Basic Installation

In the root folder `timeEvolver``, you can build the TimeEvolver with no customization using:
```
mkdir build; cd build     # create and use a build directroy
cmake ..                  # configuration reading the Cmake script
cmake --build .           # compilation and linking (or type "make")
```
This will create three folder in the folder ``build``:
* Example (example demonstrating the application of TimeEvolver to a concrete system)
* Helper (library for the creation of a matrix representation)
* TimeEvolver (core functionality: libary for the timeevolution)

Note that the generated ``Makefile`` can compile these three targets independently.


# Usage 1: Examplary Program

A first option to use the program relies on a concrete example. In order to execute the corresponding program, navigate to ``cd build/Example`` and type:
```
./main
```
A set of standard values for the parameters will be used. For a list of available parameters type:
```
./main --help
```
The result of time evolution will be stored in a HDF5-file. (For the standard choice of parameters, it has the name ``ResultBlackHole_N20_Nm2_K4_C1_DeltaN2_C00.1_Cm0.1_maxT100_tol1e-07_samplingStep0.1_m40_Sim1.h5``.) It contains the expectation values of the occupation numbers of each of the modes at different times.

# Usage 2: Apply TimeEvolver to own Hamiltonian matrix

A second option to use the program arises if the user alredy has at their disposal a Hamiltonian matrix. In this case, only the classes contained in the folder  ``TimeEvolver`` are needed. The core functionality of the TimeEvolver is encapsulated in the class ``krylovTimeEvolver`` declared in the header file ``krylovTimeEvolver.h``. Its constructor has following form
```
 krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int mm, smatrix** observables, int nbObservables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm);
```
with 
* ``t`` The time interval over which the state should be time evolved
* ``Hsize`` The size of the full Hilbert space
* ``v`` The initial state that should be time evolved
* ``samplingStep`` The time interval after which the values of the observables should be determined
* ``tol`` The maximal admissible error (norm difference between result of numerical and true time evolution
* ``m`` The size of the Krylov subspaces
* ``observables`` The matrix representations of the observables that are to be sampled
* ``nbObservables`` The number of observables
* ``Ham`` The matrix representation of the Hamiltonian
* ``expFactor`` The scalar factor multiplying the Hamiltonian in the time evolution (usually -i)
* ``checkNorm`` Whether or not it should be check that the time evolved state has unit norm

The timeevolution is started with the call of the member function
```
krylovReturn* timeEvolve();
```

## Internal data formats

Besides fundamental data types the TimeEvolver defines following two additional matrix types in the header file ``matrixDataTypes.h``:
* A dense matrix ``matrix(size_t nn, size_t mm)`` with ``nn`` the number of rows and ``mm`` the number of columns 
* A sparse matrix in coordinate format ``smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, int nn, int mm)`` with ``nn`` the number of rows, ``mm`` the number of columns, ``nbV`` the number of non-zero entries, ``row`` the 
row indices of the non-zero values, ``col`` the column indices of the non-zero values, ``val`` the non-zero values of the sparse matrix  
* The structure ``krylovReturn`` for the output of the TimeEvolver. It contains the following output data types:
    * ``matrix* sampling`` A dense matrix containing either a set of expectation values corresponding to the input observables or full quantum states at the requested time intervals
    *  ``std::complex<double>* evolvedState`` The full quantum state evolved until time ``t``
    * ``double err`` An upper error bound on the numerical error
    * ``size_t n_steps`` The number of Kyrlov steps needed for the timeevolution
    * ``size_t dim`` The dimension of the Hilbert space
    * ``size_t nSamples`` The number of samples taken

## Krylov status codes
* 0: Success
* 1: Lucky breakdown
* 2: Analytic error was smaller than the estimated total round-off error. 
* 10: Error: Numerical error estimate was larger than analytic error in at least one substep. Requested accuacy was probably not achived.  
* 11: Error: Numerical error estimate was larger than analytic error. Requested accuacy was probably not achived.
* 20: Error: Numerical integration of the error integral failed.
* 100: Error: Multible error occured.

# Build Documentation

To build a local version of the documentation you have to install the program ``doxygen``:
```
sudo apt-get install doxygen graphviz
```
Then a documentation, which lists all classes and member functions, can be built by executing in the main directory
```
doxygen
```
The result can subsequently be found by opening ``doc/html/intex.html``.


