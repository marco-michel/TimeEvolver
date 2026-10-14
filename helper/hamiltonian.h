#pragma once

#include <vector>
#include <string>
#include <complex>
#include <memory>
#include <map>

#include "matrixDataTypes.h"
#include "Basis.h"


using namespace TE;

/**
 * Represents one term in the Hamiltonian
 */ 
struct opTerm 
{
	std::complex<double> coef; //Complex coefficient of the term
	std::vector<int> operations; //Types of the involved operators (creation or annihilation)
	std::vector<int> mode; //Indices of the involved modes
	std::vector<int> capacity; //Capacities of the involved modes

	/**
 	* Turns the term into an empty one
 	*/ 
	void reset()
	{
		coef = std::complex<double>(0, 0);
		operations.clear();
		mode.clear();
		capacity.clear();
	}
};

/**
 * Basic class for all Hamiltonians that can be defined in terms of creation and annihilation operators
 */ 
class Hamiltonian
{
public:
	std::vector<opTerm> hamiltonOperator;

	Hamiltonian();

	std::string toString();
	/**
	* @param basis Basis to build the matrix in
	* @param useHermitianStorage Store only the upper triangle. The Hamiltonian
	* is Hermitian, so the lower triangle carries no information of its own, and
	* leaving it out halves both the memory and the traffic of the sparse matrix
	* vector product. Only pass true for a Hamiltonian that really is Hermitian.
	*/
	std::unique_ptr<smatrix> createHamiltonMatrix(basicBasis *basis, bool useHermitianStorage = false);
	std::vector<std::unique_ptr<smatrix>> createNumberOperatorObservables(basicBasis *basis);

	std::vector<opTerm> createKineticTerms(std::vector<int> modes, double gap);

    opTerm linInteraction(int a, int b, int ca, int cb, bool conjugate, double cof);
	opTerm createNumberOperator(int mode, double cof);
	opTerm createFourPointA(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof); 
	opTerm createFourPointB(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof);
	opTerm createFourPointC(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof);
	opTerm createSixPoint(int a, int b, int c, int d, int e, int f, int ca, int cb, int cc, int cd, int ce, int cf, double cof);


protected:
	std::unique_ptr<smatrix> createMatrix(std::vector<opTerm>& op, basicBasis* basis, bool upperTriangleOnly = false);


private:
	void creationOperator(basisState *in, unsigned int mode, int c);
	void annihilationOperator(basisState *in, unsigned int mode);
	void scalarMultiplication(basisState *in, std::complex<double> scalar);


};

