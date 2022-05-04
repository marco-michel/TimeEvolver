#include <iostream>
#include <vector>
#include <fstream>

#include <matrixDataTypes.h>
#include <hamiltonian.h>
#include <krylovTimeEvolver.h>

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
    double tol = 1.0e-8; int m = 40;

    //Number of observables = number of modes
    int nbObservables = K;

    //Creating particle number conserving basis consiting of two different quantum modes.
    std::cout << "Creating basis..." << std::endl;
    basis basis(N0, K, 0, 0);

    std::cout << "Hilbert space Dimension: " << basis.numberElements << std::endl;

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
    smatrix* hamMatrix;
    hamiltonian.createHamiltonMatrix(hamMatrix, &basis);
    //Create matrices for number Operators
    std::cout << "Creating observables..." << std::endl;
    smatrix** observables;
    hamiltonian.createObservables(observables, &basis);

    //Create initial state with all particles in the first mode (N0, 0)
    basisVector init(K); init.e[0] = N0;
    //Create initial state vector
    std::complex<double>* vec = new std::complex<double>[basis.numberElements];
	//find init state in hash table
    int entry = basis.hashTable.find(init)->second;
    vec[entry].real(1.0);

    //Start of time evolution   
    std::cout << "Starting time evolution..." << std::endl;
    std::complex<double> imaginaryMinus;
    imaginaryMinus.imag(-1);

    krylovTimeEvolver timeEvolver(maxT, basis.numberElements, vec, samplingStep, tol, m, observables, nbObservables, hamMatrix, imaginaryMinus);
    krylovReturn* results = timeEvolver.timeEvolve();

    std::cout << "Finished time evolution..." << std::endl;

    //Sort observables for easier access since observable data for different observables is stored consecutively for each time step in results matrix 
    //The storage layout of krylovReturn has folowing form: obs1 (t0), obs2(t0), obs1(t1), obs2(t1), ....
    double** nicelySorted = new double*[nbObservables];
    for(int i = 0; i < nbObservables; ++i)
        nicelySorted[i] = new double[results->nSamples]; 
    
    for (int j = 0; j < nbObservables; j++)
    {
        for (unsigned int i = 0; i < results->nSamples; i++)
            nicelySorted[j][i] = (results->sampling->values + nbObservables * i + j)->real();
    }

    //Write date to csv file
    std::cout << "Writing data to file..." << std::endl;

    for(int j = 0; j != nbObservables; j++)
    {
        std::ofstream outputfile;
        outputfile.open("output"+std::to_string(j)+".csv");
        for(int i = 0; i != results->nSamples; i++)
            outputfile << nicelySorted[j][i] << ", ";

        outputfile.close();
    }
    

    //clean up
	for (int i = 0; i < nbObservables; i++)
    {
		delete[] nicelySorted[i];
        delete observables[i];
    }
    delete[] nicelySorted; delete results; delete hamMatrix; delete[] vec; delete[] observables;

    return 0;
}