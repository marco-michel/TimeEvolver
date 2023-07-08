#pragma once


#include <vector>
#include <string>
#include <complex>

#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include "mkl_types.h"
#include "mkl.h"
#include "mkl_spblas.h"

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
	void createHamiltonMatrix(smatrix* &out, basicBasis *basis);
	void createObservables(smatrix** &out, basicBasis *basis);	

	std::vector<opTerm> createKineticTerms(std::vector<int> modes, double gap);

    opTerm linInteraction(int a, int b, int ca, int cb, bool conjugate, double cof);
	opTerm createNumberOperator(int mode, double cof);
	opTerm createFourPointA(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof); 
	opTerm createFourPointB(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof);
	opTerm createFourPointC(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof);
	opTerm createSixPoint(int a, int b, int c, int d, int e, int f, int ca, int cb, int cc, int cd, int ce, int cf, double cof);


protected:
	smatrix* createMatrix(std::vector<opTerm>& op, basicBasis* basis);


private:
	void creationOperator(basisState *in, unsigned int mode, int c);
	void annihilationOperator(basisState *in, unsigned int mode);
	void scalarMultiplication(basisState *in, std::complex<double> scalar);


};

