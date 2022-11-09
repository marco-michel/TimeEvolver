#include <vector>
#include <unordered_map>

#include "boost/container_hash/hash.hpp"
#include "boost/math/special_functions/binomial.hpp"

#include "Basis.h"

basis::~basis()
{
	delete[] basisElements;
}

/**
 * Creates a hashtable which stores the basis elements
 */
void basis::createHashTable()
{
    for (int i = 0; i != numberElements; i++)
        hashTable.insert(std::make_pair(basisElements[i], i));
}

/**
 * Creates a set of basis states. For some modes there is a fixed maximal occupation number whereas there is no such restriction for others.
 * @param N The total number of particles
 * @param K The total number of modes
 * @param nbQudits The number of modes which have a fixed maximal occupation number (K-nbQuits modes will have no limit on their occupation number)
 * @param capacity The maximal occupation number for those modes which have a limit on their occupation number
 */ 
basis::basis(int N, int K, int nbQudits, int capacity) : 
	basis(N,K,nbQudits,capacity,1) {}

/**
 * Creates a set of basis states. For some modes there is a fixed maximal occupation number whereas there is no such restriction for others.
 * @param N The total number of particles
 * @param K The total number of modes
 * @param nbQudits The number of modes which have a fixed maximal occupation number (K-nbQuits modes will have no limit on their occupation number)
 * @param capacity The maximal occupation number for those modes which have a limit on their occupation number
 * @param quant The smallest difference in occupation number in any of the modes, expect the first one (e.g. 2 means that all modes but the first one can only have even occupation numbers)
 */ 
basis::basis(int N, int K, int nbQudits, int capacity, int quant)
{

	numberModes = K;

	std::vector<basisVector> b;
	int basisSize;
	if (quant == 1)
	{
		basisSize = std::pow(capacity + 1, nbQudits)*nchoosekSmart(N + K - 1 - nbQudits, N);
	}
	else
	{
		int NEff = (int) std::ceil(N/quant);
		int capacityEff = (int) std::ceil(capacity/quant);
		basisSize = std::pow(capacityEff + 1, nbQudits)*nchoosekSmart(NEff + K - 1 - nbQudits, NEff);
		basisSize = std::min(1500000,std::abs(basisSize+1));
	}
	b.reserve(basisSize);

	basisVector test = basisVector(K);

	std::vector<basisVector> oldStates(1);
	oldStates[0] = test;
	std::vector<basisVector> newStates(basisSize);
	for (unsigned int i = 0; i != newStates.size(); i++)
		newStates[i] = test;

	basisVector newState;

	int oldStateLength = 1;
	unsigned int writeIndex = 1;
	int occuNumber = 0;
	int upperBoundary;

	basisVector dummy;


	for (int k = K; k > 1; k--) 
	{
		for (int l = 0; l <= oldStateLength-1; l++) 
		{
			occuNumber = oldStates[l].occuNumber();
			if (k > K - nbQudits)
			{
				upperBoundary = std::min(N - occuNumber, capacity);
			}
			else
			{
				upperBoundary = N - occuNumber;
			}

			for (int occ = quant; occ <= upperBoundary; occ=occ+quant)
			{
				newState = oldStates[l];
				newState.e[k-1] = occ;
				if (writeIndex >= newStates.size())
					newStates.push_back(dummy);
				newStates[writeIndex] = newState;
				writeIndex++;
			}
		}
		std::vector<basisVector>::iterator iter = newStates.begin();
		oldStateLength = writeIndex;
		oldStates.assign(iter,iter+oldStateLength);
	}

	std::vector<basisVector> preBasis(oldStateLength);
	int missing;

	writeIndex = 0; 
	
	for (int l = 0; l <= oldStateLength-1; l++) 
	{
		missing = N - oldStates[l].occuNumber();
		if (K != nbQudits || missing <= capacity)
		{
			oldStates[l].e[0] = missing;
			preBasis[writeIndex] = oldStates[l];
			writeIndex++;
		}
	}

	basisElements = new basisVector[writeIndex];

	for (unsigned int i = 0; i != writeIndex; i++)
	{
		basisElements[i] = preBasis[i];
	}

	numberElements = writeIndex;
	qudits = capacity;

    createHashTable();
	   	  
};

/**
 * Creates a tensor product of two sets of basis states. The modes in the second set have a fixed maximal occupation number whereas there is no such limit for the modes in the first set.
 * @param N1 The total number of particles in the first set
 * @param K1 The total number of modes in the first set
  * @param N2 The total number of particles in the second set
 * @param K2 The total number of modes in the second set
 * @param C2 The maximal occupation number for the modes in the second set
 */
tensorBasis::tensorBasis(int N1, int K1, int N2, int K2, int C2)
{

	numberModes = K1 + K2;
	qudits = C2;

	basis* b1 = new basis(N1, K1, 0, 0);
	basis* b2 = new basis(N2, K2, K2, C2);

	int length1 = b1->numberElements;
	int length2 = b2->numberElements;

	numberElements = length1 * length2;

	if (length2 == 0)
	{
		basisElements = new basisVector[length1];
		for (int i = 0; i != length1; i++)
			basisElements[i] = b1->basisElements[i];
		return;
	}

	basisElements = new basisVector[length1*length2];
	int ind = 0;
	basisVector tmpVec(K1 + K2);
	int* tmp = new int[K1 + K2];


	for (int i = 0; i != length1; i++)
	{
		for (int j = 0; j != length2; j++)
		{
			ind = i * length2 + j;
			for (int ii = 0; ii != K1; ii++)
				tmp[ii] = b1->basisElements[i].e[ii];
			for (int jj = K1; jj != K1+K2; jj++)
				tmp[jj] = b2->basisElements[j].e[jj-K1];
			basisElements[ind].zeros(K1 + K2);
			basisElements[ind].initInt(tmp, K1 + K2);
		}
	}


	createHashTable();


	delete b1;
	delete b2;
	delete[] tmp;
}

/**
 * Creates a hashtable which stores the basis elements
 */
void tensorBasis::createHashTable()
{
	for (int i = 0; i != numberElements; i++)
		hashTable.insert(std::make_pair(basisElements[i], i));
}


