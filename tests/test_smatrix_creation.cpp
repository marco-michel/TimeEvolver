#include <gtest/gtest.h>
#include <complex>
#include <cstddef>

#include "matrixDataTypes.h"  
#include "hamiltonian.h"
#include "exampleHamiltonian.h"


// Both tests compare a freshly created matrix against a reference stored as an
// HDF5 file, so they only exist in a build that has HDF5 available.
#ifdef USE_HDF

TEST(smatrixCreation, simpleExampleMatrix)
{
    int N0 = 200; int K = 2;
    double E1 = 1; double E2 = 2; double lambda = 1;
    double maxT = 10; double samplingStep = 0.01;
    double tol = 1.0e-6; int m = 40;
    int nbObservables = K;
    basis basis(N0, K, 0, 0);

    Hamiltonian hamiltonian;
    std::vector<opTerm> HamiltonianOperator; 
    HamiltonianOperator.push_back(hamiltonian.createNumberOperator(0, E1)); 
    HamiltonianOperator.push_back(hamiltonian.createNumberOperator(1, E2)); // adds E2 * b^dagger b
    HamiltonianOperator.push_back(hamiltonian.linInteraction(0, 1, 0, 0, false, lambda)); // adds  lambda * a^dagger b 
    HamiltonianOperator.push_back(hamiltonian.linInteraction(0, 1, 0, 0, true, lambda)); // adds  lambda * b^dagger a
    
    hamiltonian.hamiltonOperator = HamiltonianOperator;
    std::unique_ptr<smatrix> hamMatrix = hamiltonian.createHamiltonMatrix(&basis);

    smatrix storedHamMatrix;
    storedHamMatrix.loadHDF5("../../output/simpleExampleHamMatrix.h5");

    EXPECT_TRUE(hamMatrix->approxEqual(storedHamMatrix));
}


TEST(smatrixCreation, blackholeMatrix)
{
    tensorBasis basis(20, 2, 2, 2*4, 1);
    exampleHamiltonian ham = exampleHamiltonian(20, 2, 12, 4, 4, 1, 1.0, 1, 1, 1.0, 1.0);
    ham.createSimplifiedHamiltonian();
    std::unique_ptr<smatrix> hamMatrix = ham.createHamiltonMatrix(&basis);

    smatrix storedHamMatrix;
    storedHamMatrix.loadHDF5("../../output/blackHoleHamiltonian.h5");

    EXPECT_TRUE(hamMatrix->approxEqual(storedHamMatrix));
}

#endif

