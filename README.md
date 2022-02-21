![example workflow](https://github.com/marco-michel/TimeEvolver/actions/workflows/cmake.yml/badge.svg)

# Introduction


This is the TimeEvolver software package. 

# Project Structure

The TimeEvolver distribution includes the following files and directories:

<pre>
README                          this file
LICENSE                         the MIT licence
cmake                           cmake files for downloading dependencies
core                            the methods for time evolution using Krylov subspace techniques
example                         concrete example to demonstrate usage of the program
helper                          methods to create the set of basis states and to compute the Hamiltonian matrix
CMakeLists                      file needed for cmake
doxyfile                        settings file for creating documentation
</pre>


# External Requirements

This software packages relies on following external libraries:
* Intel Math Kernel Library (MKL) (required)
* BOOST (required)
* HDF5 (optional)

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
sudo apt-get install intel-mkl-full libhdf5-dev libboost-program-options-dev 
```
Note: If you install Intel MKL via the package manager you will be asked if you want to make MKL your default BLAS/LAPACK library. That is not necessary, so you can choose the default answer "No". Additionally, older versions of ``cmake`` might not find the oneapi version of mkl. We therefore recommend cmake version 3.15 or newer. 

**In order for ``TimeEvolver`` to work properly, Boost version 1.75 and newer is required. Currently, however, most distributions ship packages of the boost library of version 1.74 and older. If your machine does not have at least version 1.75, we offer several solutions (including an automatic download and compilation of the Boost-library), which are described in section "Installation"**.

## Mac

The MKL library is not supported on ARM processors used in newer Apple products. Therefore, ``TimeEvolver`` currently only works on older x86 based machines.

## Windows 

There are several ways to compile ``TimeEvolver`` on Windows. The easiest option is to use WSL (Windows Subsystem for Linux), where you can install a virtual Ubuntu machine and then follow the instructions described above. 

A second option is to download binary libraries with the help of the ``vcpk`` package manager. BOOST and HDF5 can be obtained for Windows 64bit OS via
```
.\vcpkg.exe install hdf5:x64-windows
.\vcpkg.exe install boost:x64-windows
```
You need to change the library paths in the cmake file accordingly.

A third option would be to compile Boost (and HDF5) from source. Please follow the respective instructions for each library. 

# Installation

## Basic setup

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

**Note that the above procedure will fail if you do not have at least version 1.75 of Boost. Below we describe how to solve this issue.**

We remark that the generated ``Makefile`` can compile these three targets independently.

If the dependencies have been installed locally and are not accessible system-wide one also needs to set following cmake variables with the paths of the respective libraries: ``BOOST_ROOT`` ``MKL_ROOT`` ``HDF5_DIR``. 

## Special cases

### Boost version is too old

For ``TimeEvolver`` to work properly, a Boost version of 1.75 or newer is required. If you do not have it, there are two possible solutions.

The first one is an automated download and compilation of Boost, which we integrated into cmake. To use it, you only have to set the cmake varibale ``DOWNLOAD_BOOST`` to ``ON``:
```
mkdir build; cd build        
cmake -DDOWNLOAD_BOOST=ON ..                 
cmake --build .              
```
Running the cmake script with an inadequate system-wide installed boost library before executing it with the ``DOWNLOAD_BOOST`` option might set internal variables wrongly. Deleting ``CMakeCache.txt`` and executing ``cmake -DDOWNLOAD_BOOST=ON ..`` again might solve the problem. 

The second option is to manually download the Boost-library and then to compile it from source. For a system-wide installation, one can use
```
cd path/to/boost_1_76_0
./bootstrap.sh
./b2
sudo ./b2 install
```
In case a local installation you need to change the corresponding line to ``./bootstrap.sh --prefix=INSTALL_PATH`` with ``INSTALL_PATH`` being the chosen install directory for the boost library. When using a local installation, please set the cmake variable ``-DBOOST_ROOT=INSTALL_PATH`` when building ``TimeEvolver``. 

### Standalone Version of Intel MKL

If you use the standalone version of the MKL library, which is now part of the oneAPI framework, please make sure that you follow all installation steps including setting environment variables. For root installation this can be done with:
```
source /opt/intel/oneapi/setvars.sh intel64
```
Note that the variables are only set for the context of your session. For a permanent solution please visit the Intel helppage. 

# Usage

## Usage 1: Examplary Program

A first option to use the program relies on a concrete example. In order to execute the corresponding program, navigate to ``cd build/Example`` and type:
```
./main
```
A set of standard values for the parameters will be used. For a list of available parameters type:
```
./main --help
```
The result of time evolution will be stored in a HDF5-file. (For the standard choice of parameters, it has the name ``ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.1_m40_fastIntegration0.h5``.) It contains the expectation values of the occupation numbers of each of the modes at different times.

## Usage 2: Apply TimeEvolver to own Hamiltonian matrix

A second option to use the program arises if the user alredy has at their disposal a Hamiltonian matrix. In this case, only the classes contained in the folder  ``TimeEvolver`` are needed. The core functionality of the TimeEvolver is encapsulated in the class ``krylovTimeEvolver`` declared in the header file ``krylovTimeEvolver.h``. Its constructor has following form
```
 krylovTimeEvolver(double t, size_t Hsize, std::complex<double>* v, double samplingStep, double tol, int mm, smatrix** observables, int nbObservables, smatrix* Ham, std::complex<double> expFactor, bool checkNorm, bool fastIntegration)
```
with 
* ``t`` The time interval over which the state should be time evolved
* ``Hsize`` The size of the full Hilbert space
* ``v`` The initial state that should be time evolved
* ``samplingStep`` The time interval after which the values of the observables should be determined
* ``tol`` The maximal admissible error (norm difference between result of numerical and true time evolution)
* ``m`` The size of the Krylov subspaces
* ``observables`` The matrix representations of the observables that are to be sampled
* ``nbObservables`` The number of observables
* ``Ham`` The matrix representation of the Hamiltonian
* ``expFactor`` The scalar factor multiplying the Hamiltonian in the time evolution (usually -i)
* ``checkNorm`` Whether or not it should be check that the time evolved state has unit norm
* ``fastIntegration`` Whether or not the faster not-adaptive Gauss integration scheme should be used. Note that this option can not guarantee that the numerical integration was performed with the requested accuracy.

The time evolution is started with the call of the member function
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
    * ``int statusCode``Status code indicating success or (potential) failure of the numerical time evolution 

## Krylov status codes
* ``0`` Success
* ``1`` Lucky breakdown
* ``2`` Analytic error was smaller than the estimated total round-off error. However, the estimate for the round-off error was below the requested tolerance, so the result is probably trustable. 
* ``10`` Critical Warning: Numerical error estimate was larger than analytic error in at least one substep. Requested accuacy was probably not achived.  
* ``11`` Critical Warning: Total numerical error estimate was larger than the total analytic error. Requested accuacy was probably not achived.
* ``20`` Critical Warning: Numerical integration of the error integral failed.
* ``30`` Critical Warning: Norm of state vector deviates more than the given accuaracy from 1. 
* ``100`` Critical Warning: Multible errors occured.

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

	
