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
        ("Nm", po::value<int>(&Nm)->default_value(2), "Number of particles in critical sector")
        ("K", po::value<int>(&K)->default_value(4), "Number of modes in critical sector")
        ("C0", po::value<double>(&C0)->default_value(1.0), "Coupling in control sector")
        ("Cm", po::value<double>(&Cm)->default_value(1.0), "Coupling in critical sector")
        ("maxT", po::value<double>(&maxT)->default_value(10), "Simulation-time")
	    ("samplingStep", po::value<double>(&samplingStep)->default_value(0.01), "Time interval of sampling")
        ("tol", po::value<double>(&tol)->default_value(1.0e-8), "Numerical tolerance")
        ("m", po::value<int>(&m)->default_value(40), "Dimension of Krylov-Space")
        ("DeltaN", po::value<double>(&DeltaN)->default_value(12), "Distance between critical sectors")
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
    std::unique_ptr<smatrix> hamMatrix = ham.createHamiltonMatrix(&basis);
    
   
    //Create matrices for observables
    std::cout << "Creating observables..." << std::endl;
    std::vector<std::unique_ptr<smatrix>> observables = ham.createNumberOperatorObservables(&basis);
   

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

    krylovTimeEvolver timeEvolver(maxT, vec, samplingStep, tol, m, std::move(observableList), std::move(hamMatrix), alpha, fastIntegration, true);

    timeEvolver.changeLogLevel(krylovLogger::loggingLevel::DEBUG);

    krylovReturn* results = timeEvolver.timeEvolve();

    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin);

    std::cout << "Time: " << elapsed.count()*1e-9 << " sec" << std::endl;

    std::cout << "------------------------------------------------------\n" << std::endl;
    std::cout << "Number of steps: " << results->n_steps << "    error: " << results->err << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    //End of actual time evolution 
    

    std::cout << "Writing results to file..." << std::endl;
    parameter_list parameters;

    parameters.push_back(paraPush("N", N0));
    parameters.push_back(paraPush("Nm", Nm));
    parameters.push_back(paraPush("K", K));
    parameters.push_back(paraPush("C", capacity));
    parameters.push_back(paraPush("DeltaN", DeltaN));
    parameters.push_back(paraPush("C0", C0));
    parameters.push_back(paraPush("Cm", Cm));
    parameters.push_back(paraPush("maxT", maxT));
    parameters.push_back(paraPush("tol", tol));
    parameters.push_back(paraPush("samplingStep", samplingStep));
    parameters.push_back(paraPush("m", m));
    parameters.push_back(paraPush("fastIntegration", fastIntegration));

    observableList = std::move(results->observableList);

    krylovBasicObservable::saveResult(observableList, parameters, "ResultBlackHole");

    delete results; delete[] vec;

    return 0;
}
