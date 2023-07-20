#include <iostream>
#include <vector>
#include <fstream>

#include "matrixDataTypes.h"
#include "hamiltonian.h"
#include "krylovTimeEvolver.h"
#include "krylovObservables.h"

/*Example to simulate two coupled oscillators with Hamiltonian:
 H = E1 * a^dagger a + E2 * b^dagger b + lambda * (a^dagger b + b^dagger a)
 with E1, E2 the respective gaps of the first and second quantum mode a & b
 and a coupling constant lambda
 Observables are the expectation value of the number operator in each mode
*/

int main()
{
    //Model parameters: N0-total particle number, K-number of modes
    //Energy gaps E1&E2, coulings constant lambda
    int N0 = 200; int K = 2;
    double E1 = 1; double E2 = 2; double lambda = 1;

    //Simulation parameters: Evolution time: maxT, stepsize for computing observables: samplingStep,
    //tolerance: tol, Krylovspace dimension: m (around 40 is good for most cases)
    double maxT = 10; double samplingStep = 0.01;
    double tol = 1.0e-6; int m = 40;

    //Number of observables = number of modes
    int nbObservables = K;

    //Creating particle number conserving basis consiting of two different quantum modes.
    std::cout << "Creating basis..." << std::endl;
    basis basis(N0, K, 0, 0);

    //Create Hamiltonian E1 * a^dagger a + E2 * b^dagger b + lambda * (a^dagger b + b^dagger a)
    Hamiltonian hamiltonian;
    std::vector<opTerm> HamiltonianOperator; //See definitions in hamiltonian.h for argument description
    HamiltonianOperator.push_back(hamiltonian.createNumberOperator(0, E1)); // adds E1 * a^dagger a
    HamiltonianOperator.push_back(hamiltonian.createNumberOperator(1, E2)); // adds E2 * b^dagger b
    HamiltonianOperator.push_back(hamiltonian.linInteraction(0, 1, 0, 0, false, lambda)); // adds  lambda * a^dagger b 
    HamiltonianOperator.push_back(hamiltonian.linInteraction(0, 1, 0, 0, true, lambda)); // adds  lambda * b^dagger a
    
    hamiltonian.hamiltonOperator = HamiltonianOperator;

    //Create matrix representation of HamiltonianOperator
    std::cout << "Creating Hamiltonian matrix..." << std::endl;
    std::unique_ptr<smatrix> hamMatrix = hamiltonian.createHamiltonMatrix(&basis);
    //Create matrices for number Operators
    std::cout << "Creating observables..." << std::endl;
    std::vector<std::unique_ptr<smatrix>> observables = hamiltonian.createNumberOperatorObservables(&basis);

    std::vector<std::unique_ptr<krylovBasicObservable>> observableList;
 
    for (int i = 0; i != basis.numberModes; i++)
    {
        observableList.push_back(std::make_unique<krylovSpMatrixObservable>("mode" + std::to_string(i), std::move(observables[i])));
    }

    //Create initial state with all particles in the first mode (N0, 0)
    basisVector init(K); init.e[0] = N0;
    //Create initial state vector
    std::complex<double>* vec = new std::complex<double>[basis.numberElements];
	//find init state in hash table
    int entry = basis.hashTable.find(init)->second;
    vec[entry].real(1.0);

    //Start of time evolution   
    std::cout << "Starting time evolution..." << std::endl;

    krylovTimeEvolver timeEvolver(maxT, vec, samplingStep, std::move(observableList),  std::move(hamMatrix));
    krylovReturn* results = timeEvolver.timeEvolve();

    std::cout << "Finished time evolution..." << std::endl;

    //Write date to csv file
    std::cout << "Writing data to file..." << std::endl;

    observableList = std::move(results->observableList);

    std::vector<std::unique_ptr<krylovBasicObservable>>::const_iterator obsIter;


    obsIter = observableList.begin();



   std::string outputFileName = "SimpleExampleOutputOccupationNumber";

   for (; obsIter != observableList.end(); obsIter++)
   {
       std::string fileNameCSV = outputFileName + (*obsIter)->retName() + ".csv";
       std::ofstream outputfile;
       outputfile.open(fileNameCSV);
       double* exptVal = (*obsIter)->retexpectationValues();
       for (int i = 0; i != (*obsIter)->numSamples - 1; i++)
           outputfile << exptVal[i] << ", ";
       outputfile << exptVal[((*obsIter)->numSamples) - 1];
       outputfile.close();
   }


delete[] vec; delete results;
    return 0;
}