[![CMakeTests](https://github.com/marco-michel/TimeEvolver/actions/workflows/Testscmake.yml/badge.svg?event=pull_request)](https://github.com/marco-michel/TimeEvolver/actions/workflows/Testscmake.yml)

# Introduction

*TimeEvolver* is an open-source software package for computational physics. Its purpose is to perfrom **"numerical simulations"**, i.e. to **compute real-time evolution in a generic quantum system.**  Relying on known Krylov subspace techniques, the software tackles the problem of calculating the quantum state at a given time $t$ by evaluating the product\
$\qquad\qquad\exp(-iHt)v$ ,\
where $H$ is the Hamiltonian (represented as large sparse and Hermitian matrix), *v* corresponds to the initial state vector and $i$ is the imaginary unit. The precision of the result is adjustable and on a standard notebook, *TimeEvolver* allows to compute time evolution in **any Hilbert space of dimension up to 10<sup>7</sup>**.

*TimeEvolver* is designed to be **easily applicable to your physics model**. In particular, the software package includes routines for sampling observables and for deriving the Hamiltonian matrix $H$ from a more abstract representation of the Hamiltonian operator. Moreover, convenient output methods, concrete examples and a documentation are provided.

*TimeEvolver* was considerably improved with the **new version 2.0**. Please note that the **API changed** as compared to versions 1.x (in particular of the core-method krylovTimeEvolver; see below), so existing projects may need to be updated. With the new release, *TimeEvolver* is available on **all operating systems**, although we still recommend Linux.

A detailed description of *TimeEvolver* can be found in:

M. Michel and S. Zell, *TimeEvolver*: A Program for Time Evolution With Improved Error Bound, [Comput. Phys. Commun. 277 (2022) 10837](https://doi.org/10.1016/j.cpc.2022.108374), [arXiv:2205.15346](https://arxiv.org/abs/2205.15346).

**If you use *TimeEvolver*, please cite the above paper.**

If you have any questions or experience technical problems, please feel free to **contact us**:
- <michelma@post.bgu.ac.il>
- <sebastian.zell@uclouvain.be>

# Project Structure

The *TimeEvolver* distribution includes the following files and directories:

<pre>
README                          this file
LICENSE                         the MIT licence
cmake                           cmake files for downloading dependencies
core                            the methods for time evolution using Krylov subspace techniques
example                         concrete examples to demonstrate usage of the program
helper                          methods to create the set of basis states and to compute the Hamiltonian matrix
output				output of the example
CMakeLists                      file needed for cmake
doxyfile                        settings file for creating documentation
</pre>


# External Requirements

This software packages relies on following external libraries:
* BLAS (required) - MKL recommended
* LAPACK (required) - usually included in BLAS framework
* BOOST (required)
* HDF5 (optional)

On top of that it uses the build managing software Cmake.

## Linux 

We recommend using Linux, or a virtual Linux machine, for *TimeEvolver*. The above libraries and software packages are usually available via a package manager like ``apt``. For example on Ubuntu can be installed via 
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

It is also possible to use another BLAS library instead of Intel MKL, e.g. OpenBLAS. However, the performance will most likely be worse due to lack of optimized sparse BLAS operations. In this case use
```
sudo apt-get install libopenblas-dev liblapacke-dev libhdf5-dev libboost-program-options-dev 
```

**In order for *TimeEvolver* to work properly, Boost version 1.75 and newer is required. Currently, however, most distributions ship packages of the boost library of version 1.74 and older. If your machine does not have at least version 1.75, we offer several solutions (including an automatic download and compilation of the Boost-library), which are described in section "Installation"**.

## Mac (experimental)
With version 2.0 *TimeEvolver* will also be available on Mac. The build-in ``Accelerate`` framework is used for ``BLAS`` calls. 
However, since we only have limited testing opportunities for this platform we still consider this feature experimental. For instance, ``cmake`` is not able to locate the include directory for the ``Accelerate`` framework reliably. The exact location seems also to change with different OS/xcode versions and needs to be passed by hand therefore in the ``CMakeLists.txt`` file in the root directory. Please change the following line according to your setup in case ``Accelerate.h`` is located somewhere else. 
```
include_directories(/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Headers)
```

## Windows 

There are several ways to compile *TimeEvolver* on Windows. The easiest option is to use WSL (Windows Subsystem for Linux), where you can install a virtual Ubuntu machine and then follow the instructions described above. We provide a short guide to install WSL below. 

A second option is to download binary libraries with the help of the ``vcpk`` package manager. BOOST and HDF5 can be obtained for Windows 64bit OS via
```
.\vcpkg.exe install hdf5:x64-windows
.\vcpkg.exe install boost:x64-windows
```
You need to change the library paths in the cmake file accordingly.

A third option would be to compile Boost (and HDF5) from source. Please follow the respective instructions for each library. 

# Compilation

There are two possible approaches to compile *TimeEvolver*. One can either build it locally or install it system-wide: 

**Note that both procedures will fail if you do not have at least version 1.75 of Boost. Further below we describe how to solve this issue.**

## Basic setup without installation (the easiest)

In the root folder ``timeEvolver``, you can build *TimeEvolver* with no customization using:
```
mkdir build; cd build     # create and use a build directroy
cmake ..                  # configuration reading the Cmake script
cmake --build .           # compilation and linking (or type "make")
```
This will create three folder in the folder ``build``:
* Example (example demonstrating the application of *TimeEvolver* to two concrete systems)
* Helper (library for the creation of a matrix representation)
* TimeEvolver (core functionality: libary for the timeevolution)

## Basic setup with installation

To install *TimeEvolver* to the path ``TIMEEVOLVER_INSTALL_PATH`` set the cmake variable ``CMAKE_INSTALL_PREFIX`` accordingly in the configuration step. If this variable remains unset a system folder will be chosen as installation path by cmake. 
```
mkdir build; cd build     					          # create and use a build directroy
cmake -DCMAKE_INSTALL_PREFIX=TIMEEVOLVER_INSTALL_PATH ..                  # configuration reading the Cmake script
cmake --build .          						  # compilation and linking (or type "make")
make install								  # installation to TIMEEVOLVER_INSTALL_PATH 
```

Now one can proceed as follows to compile a file that uses *TimeEvolver*. For example, if we want to compile ``simpleExample.cpp``, then we navigate to the folder where it is located and execute
```
 g++ simpleExample.cpp -I TIMEEVOLVER_INSTALL_PATH/TimeEvolver/include/ -I /opt/intel/oneapi/mkl/latest/include/ -I /usr/include/mkl -L TIMEEVOLVER_INSTALL_PATH/TimeEvolver/lib/ -lTimeEvolver -lHelper -Wl,-rpath=TIMEEVOLVER_INSTALL_PATHTimeEvolver/lib/ -o simpleExample
```
In the above, ``TIMEEVOLVER_INSTALL_PATH`` has to be replaced (three times) by the path where *TimeEvolver* was installed. Moreover, the include directoy for the MKL header might need to be adjusted. 


## Dependencies

If the dependencies have been installed locally and are not accessible system-wide one also needs to set following cmake variables with the paths of the respective libraries: ``BOOST_ROOT`` ``MKL_ROOT`` ``HDF5_DIR``. 

## Special cases

### Boost version is too old

For *TimeEvolver* to work properly, a Boost version of 1.75 or newer is required. If you do not have it, there are two possible solutions.

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
In case a local installation you need to change the corresponding line to ``./bootstrap.sh --prefix=INSTALL_PATH`` with ``INSTALL_PATH`` being the chosen install directory for the boost library. When using a local installation, please set the cmake variable ``-DBOOST_ROOT=INSTALL_PATH`` when building *TimeEvolver*. 

### Standalone Version of Intel MKL

If you use the standalone version of the MKL library, which is now part of the oneAPI framework, please make sure that you follow all installation steps including setting environment variables. For root installation this can be done with:
```
source /opt/intel/oneapi/setvars.sh intel64
```
Note that the variables are only set for the context of your session. For a permanent solution please visit the Intel helppage.

### Windows with WSL

In the following we provide a short guideline to install Windows Subsystems for Linux (WSL). For details please see the official WSL documentation. Note that WSL is only available on Windows 10 version 2004 and higher or Windows 11. Open a powershell command terminal and type
```
wsl --install
```
to install WSL. A restart is usually required after this step. To install for example ``Ubuntu`` type in a powershell or command prompt
```
wsl --install -d ubuntu
```
Note that there are many other distributions to choose from as well as newer versions of ubuntu. We again refer to the official documentation for more information. After installations you need to choose a user name and password. A Linux command terminal can be opened via the corresponding app Windows creates automatically (named after your distributions) or typing the name of the distribution in another terminal. 


## Testing

To test if the compilation was successful and the program works as intended please type
```
ctest
```
after compilation.

# Usage

## Usage 1: Examplary Program I

A first option to use the program relies on a concrete example. For this purpose we provide the code to analyze the model studied in

G. Dvali, L. Eisemann, M. Michel and S. Zell, *Black Hole Metamorphosis and Stabilization by Memory Burden*, [Phys. Rev. D, 102 (2020) 10, 103523](https://doi.org/10.1103/PhysRevD.102.103523), [arXiv:2006.00011](https://arxiv.org/abs/2006.00011).

In order to execute the corresponding program, navigate to ``cd build/Example`` and type:
```
./main
```
A set of standard values for the parameters will be used. For a list of available parameters type:
```
./main --help
```
The result of time evolution will be stored in a HDF5-file. For the standard choice of parameters, it has the name ``ResultBlackHole_N20_Nm2_K4_C1_DeltaN12_C01_Cm1_maxT10_tol1e-08_samplingStep0.01_m40_fastIntegration0.h5``. It contains the expectation values of the occupation numbers of each of the modes at different times. (If HDF5 is not installed, the result will instead be written in .csv-files.) Additionally, we provide the exemplary Mathematica-notebook ``analysisOutputData.nb`` to visualize the output data,
which is located in the example output folder ``output``.

## Usage 2: Examplary Program II

Furthermore, we provide a second simpler example of two coupled oscillators. In order to execute it, navigate to ``cd build/Example`` and type:
```
./SimpleExample
```
The expectation values of the occupation numbers of the two oscillators at different times will be written in the files ``SimpleExampleOutputOccupationNumber0.csv`` and ``SimpleExampleOutputOccupationNumber1.csv``.

## Usage 3: Apply TimeEvolver to own Hamiltonian matrix

A third option to use the program arises if the user already has at their disposal a Hamiltonian matrix. In this case, only the classes contained in the folder  ``core`` are needed. The main functionality of *TimeEvolver* is encapsulated in the class ``krylovTimeEvolver`` declared in the header file ``krylovTimeEvolver.h``. Its constructor has following form
```
krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, std::vector<std::unique_ptr<krylovBasicObservable>>  observables, std::unique_ptr<smatrix> Ham, double expFactor, double tol, int mm, bool fastIntegration, bool progressBar);
```
with 
* ``t`` The time interval over which the state should be time evolved
* ``v`` The initial state that should be time evolved
* ``samplingStep`` The time interval after which the values of the observables should be determined
* ``observables`` Vector of unique_ptrs of derived classes of krylovBasicObservable.
* ``Ham`` A unique_ptr to a matrix representation of the Hamiltonian
* ``expFactor`` The scalar factor multiplying the Hamiltonian in the time evolution (default: 1)
* ``tol`` The maximal admissible error (2-norm difference between result of numerical and true time evolution) (default: 1e-6)
* ``mm`` The size of the Krylov subspaces (default: 40)
* ``fastIntegration`` Whether or not the faster not-adaptive Gauss integration scheme should be used. Note that this option can not guarantee that the numerical integration was performed with the requested accuracy. (defaut: false)
* ``progressBar`` Whether or not a progressbar should be printed on the screen to indicate the progress.  (defaut: false)

Additionally, we provide a simplified constructor which sets all meta parameters to the above mentioned default values and only requires obligatory arguments corresponding to the physical system to be simulated. 
```
krylovTimeEvolver(double t, std::complex<double>* v, double samplingStep, std::vector<std::unique_ptr<krylovBasicObservable>> observables, std::unique_ptr<smatrix> Ham);
```

The time evolution is started with the call of the member function
```
krylovReturn* timeEvolve();
```

## Usage 4: Recent research project

*TimeEvolver* was used in a recent research project, all details can be found in the repository

[QuantumBreaking-TimeScales](https://github.com/marco-michel/QuantumBreaking-TimeScales)

Note the the previous version 1.4.1 was used there. (So small modifications would be needed before *TimeEvolver* 2.0 can be applied.)

## Internal data formats

Besides fundamental data types, *TimeEvolver* defines the following two additional matrix types in the header file ``matrixDataTypes.h``:
* A dense matrix ``matrix(size_t nn, size_t mm)`` with ``nn`` the number of rows and ``mm`` the number of columns 
* A sparse matrix in coordinate format ``smatrix(std::complex<double>* val, size_t* col, size_t* row, size_t nbV, int nn, int mm)`` with ``nn`` the number of rows, ``mm`` the number of columns, ``nbV`` the number of non-zero entries, ``row`` the 
row indices of the non-zero values, ``col`` the column indices of the non-zero values, ``val`` the non-zero values of the sparse matrix
* A base class ``krylovBasicObservable(const std::string& name)`` for handling observables, computing and storeing expectation values given a state vector as well as streamlining output to a file. There are several derivated classes to sepecalized to the cases where the observable is given by a vector, a sparse matrix or a dense matrix.
* The structure ``krylovReturn`` for the output of the *TimeEvolver*. It contains the following output data types:
    *  ``std::complex<double>* evolvedState`` The full quantum state evolved until time ``t``
    * ``double err`` An upper error bound on the numerical error
    * ``size_t n_steps`` The number of Kyrlov steps needed for the timeevolution
    * ``size_t dim`` The dimension of the Hilbert space
    * ``size_t krylovDim`` The dimension of the used Krylov space
    * ``int statusCode``Status code indicating success or (potential) failure of the numerical time evolution
    * ``std::vector<std::unique_ptr<krylovBasicObservable>> observableList`` Vector of observables each storing its own expectation values.
    * ``std::unique_ptr<smatrix>`` Hamiltonian matrix (note that entries might have been reordered)

## Krylov status codes
* ``0`` Success
* ``1`` Lucky breakdown
* ``2`` Analytic error was smaller than the estimated total round-off error. However, the estimate for the round-off error was below the requested tolerance, so the result is probably trustable. 
* ``3`` The TimeEvolution was aborted by request of an observable. KrylovReturn only contains results up until exception was thrown. Note that warnings that have been thrown up until this point are not retrievable. 
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
The result can subsequently be found by opening ``doc/html/index.html``.

	
