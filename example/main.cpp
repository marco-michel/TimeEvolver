#include <unordered_map>
#include <complex>
#include <iostream>
#include <string>
#include <memory>
#include <chrono>


#include <boost/program_options.hpp>
namespace po = boost::program_options;



#include "matrixDataTypes.h"
#include "krylovTimeEvolver.h"
#include "Basis.h"
#include "hamiltonian.h"
#include "exampleHamiltonian.h"
#include "krylovObservables.h"
#include "parameter.h"
#include "krylovHelper.h"

using namespace TE;

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
    double tol; int m;
    double DeltaN; int capacity; 
    bool fastIntegration;
    
    po::options_description desc("Allowed options");
    po::variables_map vm;
    
    try {
        
        desc.add_options()
        ("help", "produce help message")
        ("N0", po::value<int>(&N0)->default_value(20), "Number of particles in control sector")
        ("Nm", po::value<int>(&Nm)->default_value(3), "Number of particles in critical sector")
        ("K", po::value<int>(&K)->default_value(6), "Number of modes in critical sector")
        ("C0", po::value<double>(&C0)->default_value(0.0003), "Coupling in control sector")
        ("Cm", po::value<double>(&Cm)->default_value(0.065), "Coupling in critical sector")
        ("maxT", po::value<double>(&maxT)->default_value(500), "Simulation-time")
	    ("samplingStep", po::value<double>(&samplingStep)->default_value(20), "Time interval of sampling")
        ("tol", po::value<double>(&tol)->default_value(1.0e-6), "Numerical tolerance")
        ("m", po::value<int>(&m)->default_value(40), "Dimension of Krylov-Space")
        ("DeltaN", po::value<double>(&DeltaN)->default_value(12), "Distance between critical sectors")
        ("capacity", po::value<int>(&capacity)->default_value(1), "Capacity of cirtial modes")
        ("fastIntegration", po::value<bool>(&fastIntegration)->default_value(true), "Use faster and less accurate integration")
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

    mkl_set_num_threads(2);

    std::cout << "TimeEvolver Example" << std::endl;
    
    //Create basis
    tensorBasis basis(N0, 2, Nm, 2*K, capacity);
    std::cout << "Created basis with " << basis.numberElements << " elements ..." << std::endl;
    
    //Create Hamiltonian
    exampleHamiltonian ham = exampleHamiltonian(N0, Nm, DeltaN, K, K, capacity, C0, 1, 1, Cm, Cm);
    ham.createSimplifiedHamiltonian();
    //std::cout << ham.toString() << std::endl;
    
   //Create Hamiltonian matrix
    std::unique_ptr<smatrix> hamMatrix = ham.createHamiltonMatrix(&basis);
    std::cout << "Created Hamiltonian matrix..." << std::endl;
    
   
    //Create matrices for observables
    std::vector<std::unique_ptr<smatrix>> observables = ham.createNumberOperatorObservables(&basis);
    std::cout << "Created observables..." << std::endl;
   

    //Create initial state
    basisVector init = ham.createInitState();
    std::complex<double>* vec = new std::complex<double>[basis.numberElements];
	//find init state in hash table
    int entry = basis.hashTable.find(init)->second;
    vec[entry].real(1.0);
   
    //Start of actual time evolution   
    std::cout << "Starting time evolution..." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<std::unique_ptr<krylovBasicObservable>> observableList;
 
    for (int i = 0; i != basis.numberModes; i++)
    {
        observableList.push_back(std::make_unique<krylovSpMatrixObservable>("mode" + std::to_string(i), std::move(observables[i])));
    }
  
    double alpha = 1.0;

    krylovTimeEvolver timeEvolver(maxT, vec, samplingStep, std::move(observableList), std::move(hamMatrix), alpha, tol, m, fastIntegration, true);

    timeEvolver.changeLogLevel(krylovLogger::loggingLevel::DEBUG);

    krylovReturn* results = timeEvolver.timeEvolve();

    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin);

    std::cout << "Time: " << elapsed.count()*1e-9 << " sec" << std::endl;

    std::cout << "------------------------------------------------------\n" << std::endl;
    std::cout << "Number of steps: " << results->n_steps << "    error: " << results->err << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    //End of actual time evolution 
    

    parameter_list parameters;

    parameters.push_back(paraPush("N", true, N0));
    parameters.push_back(paraPush("Nm", true, Nm));
    parameters.push_back(paraPush("K", true, K));
    parameters.push_back(paraPush("C", true, capacity));
    parameters.push_back(paraPush("DeltaN", true, DeltaN));
    parameters.push_back(paraPush("C0", true, C0));
    parameters.push_back(paraPush("Cm", true, Cm));
    parameters.push_back(paraPush("maxT", true, maxT));
    parameters.push_back(paraPush("tol", true, tol));
    parameters.push_back(paraPush("samplingStep", true, samplingStep));
    parameters.push_back(paraPush("m", true, m));
    parameters.push_back(paraPush("fastIntegration", true, fastIntegration));
    parameters.push_back(paraPush("TE_MAJOR_VERSION", false, TIMEEVOLVER_VERSION / 100));
    parameters.push_back(paraPush("TE_MINOR_VERSION", false, TIMEEVOLVER_VERSION % 100));

    observableList = std::move(results->observableList);

    saveResult(observableList, parameters, "ResultBlackHole");
    std::cout << "Results have been saved to file." << std::endl;

    delete results; delete[] vec;

    return 0;
}
