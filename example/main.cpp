#include <unordered_map>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include <chrono>

#undef I
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <boost/program_options.hpp>
#include <H5Cpp.h>
#include <mkl_types.h>
#include <mkl.h>

#include "matrixDataTypes.h"
#include "krylovTimeEvolver.h"
#include "Basis.h"
#include "hamiltonian.h"
#include "exampleHamiltonian.h"


namespace po = boost::program_options;
using namespace H5;


//INPUT:
//    int N0; int Nm; int K;
//  double C0; double Cm
//  double maxT; double samplingStep; 
//   double tol; int m; int numThreads; int DeltaN; int capacity;  bool fastIntegration
//  in total, argc can take up to 13 arguments
int main(int argc, char* argv[])
{
    
    //Determine parameters
    int N0; int Nm;
    int K; 
    double C0; double Cm;
    double maxT; double samplingStep;
    double tol; int m; int numThreads;
     double DeltaN; int capacity; 
     bool fastIntegration;
    
    po::options_description desc("Allowed options");
    po::variables_map vm;
    
    try {
        
        desc.add_options()
        ("help", "produce help message")
        ("N0", po::value<int>(&N0)->default_value(20), "Number of particles in control sector")
        ("Nm", po::value<int>(&Nm)->default_value(2), "Number of particles in critical sector")
        ("K", po::value<int>(&K)->default_value(4), "Number of modes in critical sector")
        ("C0", po::value<double>(&C0)->default_value(0.1), "Coupling in control sector")
        ("Cm", po::value<double>(&Cm)->default_value(0.1), "Coupling in critical sector")
        ("maxT", po::value<double>(&maxT)->default_value(100), "Simulation-time")
	    ("samplingStep", po::value<double>(&samplingStep)->default_value(0.1), "Time interval of sampling")
        ("tol", po::value<double>(&tol)->default_value(1.0e-7), "Numerical tolerance")
        ("m", po::value<int>(&m)->default_value(40), "Dimension of Krylov-Space")
        ("threads", po::value<int>(&numThreads)->default_value(2), "Number of OpenMP Threads for Intel MKL")
        ("DeltaN", po::value<double>(&DeltaN)->default_value(2), "Distance between critical sectors")
        ("capacity", po::value<int>(&capacity)->default_value(1), "Capacity of cirtial modes")
        ("fastIntegration", po::value<bool>(&fastIntegration)->default_value(false), "Use faster and less accurate integration")
        ;
    
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
    } catch (const boost::program_options::required_option & e) {
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        } else {
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }
    
    
    int nbObservables = 2*K + 2;

    mkl_set_num_threads(numThreads);
    //end determine parameters


    std::cout << "TimeEvolver Example" << std::endl;
    
    //Create basis
    std::cout << "Creating basis..." << std::endl;
    tensorBasis basis(N0, 2, Nm, 2*K, capacity);
    
    //Create Hamiltonian
    exampleHamiltonian ham = exampleHamiltonian(N0, Nm, DeltaN, K, K, capacity, C0, 1, 1, Cm, Cm);
    ham.createSimplifiedHamiltonian();
    //std::cout << ham.toString() << std::endl;
    
   //Create Hamiltonian matrix
    std::cout << "Creating Hamiltonian matrix..." << std::endl;
    smatrix* hamMatrix;
    ham.createHamiltonMatrix(hamMatrix, &basis);
    //Create matrices for observables
    std::cout << "Creating observables..." << std::endl;
    smatrix** observables;
    ham.createObservables(observables, &basis);

    //Create initial state
    basisVector init = ham.createInitState();
    std::complex<double>* vec = new std::complex<double>[basis.numberElements];
	//find init state in hash table
    int entry = basis.hashTable.find(init)->second;
    vec[entry].real(1.0);
   
    //Start of actual time evolution   
    std::cout << "Starting time evolution..." << std::endl;
    std::complex<double> imaginaryMinus;
    imaginaryMinus.imag(-1);

    auto begin = std::chrono::high_resolution_clock::now();

    krylovTimeEvolver timeEvolver(maxT, basis.numberElements, vec, samplingStep, tol, m, observables, nbObservables, hamMatrix, imaginaryMinus, true, fastIntegration);
    krylovReturn* results = timeEvolver.timeEvolve();
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin);

    std::cout << "Time: " << elapsed.count()*1e-9 << " sec" << std::endl;

    std::cout << "------------------------------------------------------\n" << std::endl;
    std::cout << "Number of steps: " << results->n_steps << "    error: " << results->err << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    //End of actual time evolution    
    
    //Create string with information about input parameters
        std::cout << "Writing results to file..." << std::endl;
       const int nbOutputParameters = 12;
    std::ostringstream outputnumbers[nbOutputParameters];
    outputnumbers[0] << std::fixed << std::setprecision(0) << N0;
    outputnumbers[1] << std::fixed << std::setprecision(0) << Nm;
    outputnumbers[2] << std::fixed << std::setprecision(0) <<K;
    outputnumbers[3] << std::fixed << std::setprecision(0) << capacity;
    outputnumbers[4] << DeltaN;
    outputnumbers[5] << C0; 
    outputnumbers[6] << Cm;
    outputnumbers[7] << std::fixed << std::setprecision(0) << maxT;
    outputnumbers[8] << tol;
    outputnumbers[9] << samplingStep;
    outputnumbers[10] << std::fixed << std::setprecision(0) << m;
     outputnumbers[11] << fastIntegration;
    
    std::string obligatoryInfo = "_N" + outputnumbers[0].str() + "_Nm" + outputnumbers[1].str() + "_K"
    + outputnumbers[2].str() + "_C" + outputnumbers[3].str();
    
    std::string furtherInfo = "_DeltaN" + outputnumbers[4].str() + "_C0" + outputnumbers[5].str() + "_Cm" + outputnumbers[6].str()
    + "_maxT" + outputnumbers[7].str() + "_tol" + outputnumbers[8].str() + "_samplingStep" + outputnumbers[9].str() + "_m" + outputnumbers[10].str() + "_fastIntegration" + outputnumbers[11].str() + "_Sim1";
   
    std::string fileNameH5 = "ResultBlackHole" + obligatoryInfo + furtherInfo + ".h5";
    //End create string
    
    //Write observables to file
    double** nicelySorted = new double*[nbObservables];
    for(int i = 0; i < nbObservables; ++i)
        nicelySorted[i] = new double[results->nSamples]; 
    
    for (int j = 0; j < nbObservables; j++)
    {
        for (unsigned int i = 0; i < results->nSamples; i++)
        {
            nicelySorted[j][i] = (results->sampling->values + nbObservables * i + j)->real();
        }
    }

    H5File fileHh(fileNameH5.c_str(), H5F_ACC_TRUNC);

    hsize_t dims[1] = { 1 };
    
    
    for (int j = 0; j < nbObservables; j++)
    {
        int NX = results->nSamples;
        const int RANK = 1;
        hsize_t dimsf[RANK];
        dimsf[0] = NX;
        DataSpace dataspace( RANK, dimsf );
        FloatType datatype( PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );
		std::string modeNumber = "mode" + std::to_string(j);
        DataSet dataset = fileHh.createDataSet(modeNumber.c_str(), datatype,
                                              dataspace );
        dataset.write( nicelySorted[j], PredType::NATIVE_DOUBLE );

        
        //write all parameters as attributes to file
        
        DataSpace** attr_dataspace = new DataSpace*[nbOutputParameters];
        Attribute** attributes = new Attribute*[nbOutputParameters];
        
        for (int i = 0; i < nbOutputParameters; ++i)
        {
            attr_dataspace[i] = new DataSpace (1, dims );
            attributes[i] = new Attribute();
        }
        
        
        
        *attributes[0] = dataset.createAttribute( "DeltaN", PredType::NATIVE_DOUBLE, *attr_dataspace[0]); attributes[0]->write(PredType::NATIVE_DOUBLE, &DeltaN);
        *attributes[1] = dataset.createAttribute( "C0", PredType::NATIVE_DOUBLE, *attr_dataspace[1]); attributes[1]->write(PredType::NATIVE_DOUBLE, &C0);
        *attributes[2] = dataset.createAttribute( "Cm", PredType::NATIVE_DOUBLE, *attr_dataspace[2]); attributes[2]->write(PredType::NATIVE_DOUBLE, &Cm);
        *attributes[3] = dataset.createAttribute( "MaxT", PredType::NATIVE_DOUBLE, *attr_dataspace[3]); attributes[3]->write(PredType::NATIVE_DOUBLE, &maxT);
        *attributes[4] = dataset.createAttribute( "Tolerance", PredType::NATIVE_DOUBLE, *attr_dataspace[4]); attributes[4]->write(PredType::NATIVE_DOUBLE, &tol);
        *attributes[5] = dataset.createAttribute( "SamplingStep", PredType::NATIVE_DOUBLE, *attr_dataspace[5]); attributes[5]->write(PredType::NATIVE_DOUBLE, &samplingStep);
        *attributes[6] = dataset.createAttribute( "N0", PredType::NATIVE_INT, *attr_dataspace[6]); attributes[6]->write(PredType::NATIVE_INT, &N0);
        *attributes[7] = dataset.createAttribute( "Nm", PredType::NATIVE_INT, *attr_dataspace[7]); attributes[7]->write(PredType::NATIVE_INT, &Nm);
        *attributes[8] = dataset.createAttribute( "K", PredType::NATIVE_INT, *attr_dataspace[8]); attributes[8]->write(PredType::NATIVE_INT, &K);
        *attributes[9] = dataset.createAttribute( "Capacity", PredType::NATIVE_INT, *attr_dataspace[9]); attributes[9]->write(PredType::NATIVE_INT, &capacity);
        *attributes[10] = dataset.createAttribute( "m", PredType::NATIVE_INT, *attr_dataspace[10]); attributes[10]->write(PredType::NATIVE_INT, &m);
        *attributes[11] = dataset.createAttribute( "fastIntegration", PredType::NATIVE_INT, *attr_dataspace[11]); attributes[11]->write(PredType::NATIVE_INT, &fastIntegration);
        
        for (int i = 0; i < nbOutputParameters; ++i)
        {
            delete attr_dataspace[i];
            delete attributes[i];
        }


        delete[] attr_dataspace;
        delete[] attributes;
        dataset.close();

    }
    //end write obserables to file


    //clean up
	for (int i = 0; i < nbObservables; i++)
		delete[] nicelySorted[i];
	delete[] nicelySorted;
    
    delete results; 

    delete hamMatrix;
    delete[] vec;

    for (int i = 0; i != nbObservables; i++)
    {
        delete observables[i];
    }
    if(nbObservables > 0 )
        delete[] observables;
    
    return 0;
}
