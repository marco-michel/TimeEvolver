#include <complex>
#include <vector>
#include <map>

#include "matrixDataTypes.h"
#include "hamiltonian.h"
#include "Basis.h"


/**
 * Creates the kinetic terms
 * @param modes Indices of the modes for which the kinetic terms should be created
 * @param gap Energy gap of the modes
 * @return Returns operators corresponding to the kinetic terms
 */
std::vector<opTerm> Hamiltonian::createKineticTerms(std::vector<int> modes, double gap)
{
	std::complex<double> eins = std::complex<double>(gap, 0.0);
	std::vector<opTerm> ret;
	opTerm tmp;

	for (std::vector<int>::iterator iter = modes.begin(); iter != modes.end(); iter++)
	{
		tmp.coef = eins;
		tmp.capacity.push_back(0);
		tmp.capacity.push_back(0);
		tmp.mode.push_back(*iter);
		tmp.mode.push_back(*iter);
		tmp.operations.push_back(-1);
		tmp.operations.push_back(1);

		ret.push_back(tmp);
		tmp.reset();
	}
	return ret;
}

/**
 * Creates an interaction term involving one creation and one annihilation operator
 * @param a Index of the modes on which the first operator acts
 * @param b Index of the modes on which the second operator acts
 * @param ca Maximum occupation number of the modes on which the first operator acts (0 if no maximum)
 * @param cb Maximum occupation number of the modes on which the second operator acts (0 if no maximum)
 * @param conjugate If true, return the Hermitian conjugate of what would have otherwise been returned
 * @param cof Complex coefficient of the operator
 * @return Returns operator corresponding to the bilinear interaction
 */
opTerm  Hamiltonian::linInteraction(int a, int b, int ca, int cb, bool conjugate, double cof)
{
    opTerm tmp;
    
    tmp.coef.real(cof);
    tmp.operations.push_back(-1);
    tmp.operations.push_back(1);
    if(!conjugate)
    {
        tmp.capacity.push_back(cb);
        tmp.capacity.push_back(ca);
        tmp.mode.push_back(b);
        tmp.mode.push_back(a);
    }
    else
    {
        tmp.capacity.push_back(ca);
        tmp.capacity.push_back(cb);
        tmp.mode.push_back(a);
        tmp.mode.push_back(b);
    }
    return tmp;
}

/**
 * Creates a number operator
 * @param mode Index of the mode for which the number operator should be created
 * @param cof Complex coefficient of the operator
 * @return Returns number operator
 */
opTerm Hamiltonian::createNumberOperator(int mode, double cof)
{
	opTerm tmp;

	tmp.coef.real(cof);
	tmp.capacity.push_back(0);
	tmp.capacity.push_back(0);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);
	tmp.mode.push_back(mode);
	tmp.mode.push_back(mode);
	return tmp;
}

/**
 * Creates an interaction term consisting of 4 operators, namely in order creation, annihilation, creation, annihilation
 * @param a Index of the modes on which the first operator acts
 * @param b Index of the modes on which the second operator acts
 * @param c Index of the modes on which the third operator acts
 * @param d Index of the modes on which the fourth operator acts
 * @param ca Maximum occupation number of the modes on which the first operator acts (0 if no maximum)
 * @param cb Maximum occupation number of the modes on which the second operator acts (0 if no maximum)
 * @param cc Maximum occupation number of the modes on which the third operator acts (0 if no maximum)
 * @param cd Maximum occupation number of the modes on which the fourth operator acts (0 if no maximum)
 * @param cof Complex coefficient of the operator
 * @return Returns operator corresponding to the 4-point interaction
 */
opTerm Hamiltonian::createFourPointA(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof)
{
	opTerm tmp;

	tmp.coef.real(cof);
	tmp.capacity.push_back(cd);
	tmp.capacity.push_back(cc);
	tmp.capacity.push_back(cb);
	tmp.capacity.push_back(ca);
	tmp.mode.push_back(d);
	tmp.mode.push_back(c);
	tmp.mode.push_back(b);
	tmp.mode.push_back(a);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);

	return tmp;
}

/**
 * Creates an interaction term consisting of 4 operators, namely in order creation, creation, annihilation, annihilation
 * @param a Index of the modes on which the first operator acts
 * @param b Index of the modes on which the second operator acts
 * @param c Index of the modes on which the third operator acts
 * @param d Index of the modes on which the fourth operator acts
 * @param ca Maximum occupation number of the modes on which the first operator acts (0 if no maximum)
 * @param cb Maximum occupation number of the modes on which the second operator acts (0 if no maximum)
 * @param cc Maximum occupation number of the modes on which the third operator acts (0 if no maximum)
 * @param cd Maximum occupation number of the modes on which the fourth operator acts (0 if no maximum)
 * @param cof Complex coefficient of the operator
 * @return Returns operator corresponding to the 4-point interaction
 */
opTerm Hamiltonian::createFourPointB(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof)
{
	opTerm tmp;

	tmp.coef.real(cof);
	tmp.capacity.push_back(cd);
	tmp.capacity.push_back(cc);
	tmp.capacity.push_back(cb);
	tmp.capacity.push_back(ca);
	tmp.mode.push_back(d);
	tmp.mode.push_back(c);
	tmp.mode.push_back(b);
	tmp.mode.push_back(a);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);
	tmp.operations.push_back(1);

	return tmp;
}

/**
 * Creates an interaction term consisting of 4 operators, namely in order annihilation, annihilation, creation, creation
 * @param a Index of the modes on which the first operator acts
 * @param b Index of the modes on which the second operator acts
 * @param c Index of the modes on which the third operator acts
 * @param d Index of the modes on which the fourth operator acts
 * @param ca Maximum occupation number of the modes on which the first operator acts (0 if no maximum)
 * @param cb Maximum occupation number of the modes on which the second operator acts (0 if no maximum)
 * @param cc Maximum occupation number of the modes on which the third operator acts (0 if no maximum)
 * @param cd Maximum occupation number of the modes on which the fourth operator acts (0 if no maximum)
 * @param cof Complex coefficient of the operator
 * @return Returns operator corresponding to the 4-point interaction
 */
opTerm Hamiltonian::createFourPointC(int a, int b, int c, int d, int ca, int cb, int cc, int cd, double cof)
{
	opTerm tmp;

	tmp.coef.real(cof);
	tmp.capacity.push_back(cd);
	tmp.capacity.push_back(cc);
	tmp.capacity.push_back(cb);
	tmp.capacity.push_back(ca);
	tmp.mode.push_back(d);
	tmp.mode.push_back(c);
	tmp.mode.push_back(b);
	tmp.mode.push_back(a);
	tmp.operations.push_back(1);
	tmp.operations.push_back(1);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(-1);

	return tmp;
}

/**
 * Creates an interaction term consisting of 6 operators, namely in order creation, annihilation, creation, annihilation, creation, annihilation
 * @param a Index of the modes on which the first operator acts
 * @param b Index of the modes on which the second operator acts
 * @param c Index of the modes on which the third operator acts
 * @param d Index of the modes on which the fourth operator acts
 * @param e Index of the modes on which the fifth operator acts
 * @param f Index of the modes on which the sixth operator acts
 * @param ca Maximum occupation number of the modes on which the first operator acts (0 if no maximum)
 * @param cb Maximum occupation number of the modes on which the second operator acts (0 if no maximum)
 * @param cc Maximum occupation number of the modes on which the third operator acts (0 if no maximum)
 * @param cd Maximum occupation number of the modes on which the fourth operator acts (0 if no maximum)
 * @param ce Maximum occupation number of the modes on which the fifth operator acts (0 if no maximum)
 * @param cf Maximum occupation number of the modes on which the sixth operator acts (0 if no maximum)
 * @param cof Complex coefficient of the operator
 * @return Returns operator corresponding to the 6-point interaction
 */
opTerm Hamiltonian::createSixPoint(int a, int b, int c, int d, int e, int f, int ca, int cb, int cc, int cd, int ce, int cf, double cof)
{
	opTerm tmp;

	tmp.coef.real(cof);
	tmp.capacity.push_back(cf);
	tmp.capacity.push_back(ce);
	tmp.capacity.push_back(cd);
	tmp.capacity.push_back(cc);
	tmp.capacity.push_back(cb);
	tmp.capacity.push_back(ca);
	tmp.mode.push_back(f);
	tmp.mode.push_back(e);
	tmp.mode.push_back(d);
	tmp.mode.push_back(c);
	tmp.mode.push_back(b);
	tmp.mode.push_back(a);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);
	tmp.operations.push_back(-1);
	tmp.operations.push_back(1);

	return tmp;
}

/**
 * Creates an empty Hamiltonian
 */
Hamiltonian::Hamiltonian()
{

}


/**
 * Represents the Hamiltonian as string
 * @return string representation of the Hamiltonian
 */
std::string Hamiltonian::toString() 
{
	if (hamiltonOperator.empty())
	{
		return std::string();
	}

	std::string ret;
	std::vector<opTerm>::iterator iter = hamiltonOperator.begin();
	for (; iter != hamiltonOperator.end(); iter++)
	{
		ret += std::to_string(iter->coef.real());
		for (unsigned int i = 0; i != iter->capacity.size(); i++)
		{
			if (iter->capacity[i] == 0)
				ret += "a";
			else
				ret += "b";
			ret += std::to_string(iter->mode[i]);
			if (iter->operations[i] == 1)
				ret += "d";
		}
		ret += "\n+";
	}

	return ret;
}

/**
 * Creates the represenation of the Hamiltonian as a sparse matrix
 * @param out Resulting sparse matrix
 * @param basis Set of basis states
 */
void Hamiltonian::createHamiltonMatrix(smatrix * &out, basicBasis* basis)
{
	out = createMatrix(hamiltonOperator, basis);
}

/**
 * Creates an array of observables (represented as matrices). For each mode, an observable corresponding to the expectation value of the of the occupation of the mode is created.
 * @param out resulting list of observables
 * @param basis set of basis states
 */
void Hamiltonian::createObservables(smatrix ** &out, basicBasis* basis)
{
    int nbModes = basis->numberModes;

	out = new smatrix*[nbModes];

	for (int i = 0; i != 2; i++)
	{
		std::vector<opTerm> operators;
		operators.push_back(createNumberOperator(i, 1));
		smatrix* matrix;
		matrix = createMatrix(operators, basis);
		out[i] = matrix;
	}

	for (int i = 2; i != nbModes; i++)
	{
		std::vector<opTerm> operators;
		operators.push_back(createNumberOperator(i, 1));
		smatrix* matrix;
		matrix = createMatrix(operators, basis);
		out[i] = matrix;
	}


}

/**
 * Multiplies basis state by a scalar coefficient
 * @param in state to be multiplied
 * @param scalar scalar coefficient used to multiply
 */
void Hamiltonian::scalarMultiplication(basisState* in, std::complex<double> scalar)
{
	in->coef = in->coef * scalar;
}

/**
 * Creates the represenation of an operator as sparse matrix
 * @param op terms that sum up to form the total operator
 * @param basis set of basis states
 */
smatrix* Hamiltonian::createMatrix(std::vector<opTerm>& op, basicBasis * basis)
{

	std::complex<double> one(1, 0);

	basisState* tmp = nullptr;

	int matrixSize = basis->numberElements*op.size();

	std::complex<double>* values = new std::complex<double>[matrixSize];
	size_t* rowIndex = new size_t[matrixSize];
	size_t *columnIndex = new size_t[matrixSize];
	int Index = 0;
	std::vector<int> colTmp;
	std::vector<std::complex<double>> valTmp;
	std::map<int, std::complex<double>> perRow;
	std::map<int, std::complex<double>>::iterator perRowIter;
	std::map<int, std::complex<double>>::iterator pairDoesNotWork;
	std::pair<std::map<int,std::complex<double>>::iterator, bool> tester;


	for (int j = 0; j != basis->numberElements; j++) //loop over basis elements
	{

		std::vector<opTerm>::iterator iter = op.begin(); 



		perRow.clear();

		for (; iter != op.end(); iter++) //loop over terms in hamiltonian
		{
			tmp = new basisState(basis->basisElements[j], one);

			for (unsigned int i = 0; i != iter->operations.size(); i++) //loop over operations in each term
			{
				if (iter->operations[i] == 1)
				{
					creationOperator(tmp, iter->mode[i], iter->capacity[i]);
				}
				else if (iter->operations[i] == -1)
				{
					annihilationOperator(tmp, iter->mode[i]);
				}
			}

			scalarMultiplication(tmp, iter->coef);

			bool skip = false;
			int currentElement = -1;
			std::unordered_map<basisVector, int, basisVectorHasher>::iterator hashTableIterator;
			hashTableIterator = basis->hashTable.find(tmp->b);
			if (hashTableIterator != basis->hashTable.end())
				currentElement = hashTableIterator->second;
			else
			{
				skip = true;
			}
			if (!skip)
			{
				pairDoesNotWork = perRow.find(currentElement);
				if (pairDoesNotWork == perRow.end())
					perRow.insert(std::make_pair(currentElement, tmp->coef));
				else
					pairDoesNotWork->second += tmp->coef;
			}
			skip = false;
			delete tmp;

		}

		//write row into matrix
		perRowIter = perRow.begin();
		for (; perRowIter != perRow.end(); perRowIter++)
		{
			if (std::abs(perRowIter->second) < std::numeric_limits<double>::epsilon()) // no zeros 
				continue;
			rowIndex[Index] = j;
			columnIndex[Index] = perRowIter->first;
			values[Index] = perRowIter->second;
			Index++;
		}


	}

	int M, N;
	M = N = basis->numberElements;
	int nz = Index;

	smatrix* A = new smatrix(values, columnIndex, rowIndex, nz, N, M);

	delete[] rowIndex;
	delete[] columnIndex;
	delete[] values;

	return A;
}

/**
 * Applies creation operator to a basis state
 * @param in basis state to be modified by the operator
 * @param mode number of mode on which creation operator acts
 * @param c maximal occupation number of that mode (0 if no limit on occupation)
 */
void Hamiltonian::creationOperator(basisState *in, unsigned int mode, int c)
{
	if (mode >= in->b.length)
	{
		std::cerr << "creationOperator: mode does not exist in vector\n";
		exit(14);
	}

	int n = in->b.e[mode];
	if (c == 0)
	{
		in->coef = in->coef*sqrt(n + 1.0);
		in->b.e[mode] = n + 1;
	}
	else
	{
		if (n + 1 > c)
			in->coef = 0;
		else
		{
			in->coef = in->coef*sqrt(n + 1.0);
			in->b.e[mode] = n + 1;
		}
	}
}

/**
 * Applies annihilation operator to a basis state
 * @param in basis state to be modified by the operator
 * @param mode number of mode on which annihilation operator acts
 */
void Hamiltonian::annihilationOperator(basisState* in, unsigned int mode)
{
	if (mode >= in->b.length)
	{
		std::cerr << "creationOperator: mode does not exist in vector\n";
		exit(14);
	}
	int n = in->b.e[mode];
	if (n == 0)
		in->coef = 0;
	else
	{
		in->coef = in->coef*sqrt(n);
		in->b.e[mode] = n - 1;
	}
}


