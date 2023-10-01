#pragma once

#include <complex>
#include <vector>
#include <iostream>
#include <unordered_map>

#include "boost/container_hash/hash.hpp"
#include "boost/math/special_functions/binomial.hpp"


/**
 * A basis vector in the basis of number eigenstates. It is represented as a vector of length l, where l is the number of modes in the system. The kth entry of the vector corresponds to the occupation number of the mode k.
 */
struct basisVector
{
	int* e;
	size_t length;

	/**
 	* Creates a vector of length 0
	 */
	basisVector()
	{
		length = 0;
		e = nullptr;
	}

	/**
 	* Creates a vector in which all modes are unoccupied
 	 * @param le The length of the vector, i.e. the number of modes
	 */
	basisVector(unsigned int le)
	{
		if (le == 0)
		{
			e = nullptr;
			length = 0;
			return;
		}

		length = le;
		e = new int[le]();
	}
	
	/**
 	* Creates a vector in which all modes but one mode are unoccupied
 	 * @param le The length of the vector, i.e. the number of modes
 	  * @param pos The position of the mode which is occupied (with occupation number 1)
	 */
	basisVector(int le, int pos)
	{
		if (le <= pos)
		{
			std::cerr << "basisvector too small" << std::endl;
			exit(12);
		}
		if (le > 0)
		{
			length = le;
			e = new int[le]();
			e[pos] = 1;
		}
		else
		{
			length = 0;
			e = nullptr;
		}
	}
	
	/**
 	* Sets the occupation number of all modes to 1
	 */
	void ones()
	{
		for (unsigned int i = 0; i < length; i++)
		{
			e[i] = 1;
		}
	}

	/**
 	* Sets the occupation number of all modes to 0
	 */
	void zeros()
	{
		for (unsigned int i = 0; i < length; i++)
		{
			e[i] = 0;
		}
	}

	/**
 	* Changes the length of the vector and then sets all occupation number to 0
 	* @param n The new length of the vector
	 */
	void zeros(int n)
	{
		if (n > 0)
		{
			if (length != 0)
				delete[] e;
			length = n;
			e = new int[n];
		}
		zeros();
	}

	/**
 	* Changes the length of the vector and then sets all occupation numbers
 	* @param init An array which specifies the occupation numbers
 	* @param le The new length of the vector
	 */
	void initInt(int * init, unsigned int le)
	{
		if (length == le)
		{
			for (unsigned int i = 0; i != le; i++)
				e[i] = init[i];
		}
	}
		
	/**
 	* Computes the total occupation number in all modes
	 */
	int occuNumber()
	{
		int ret = 0;
		if (length == 0)
			return ret;
		for (unsigned int i = 0; i != length; i++)
			ret += e[i];
		return ret;
	}

	/**
 	* Creates a new basis vector by copying an existing one
 	* @param old_obj The existing object to be copied
	 */
	basisVector(const basisVector &old_obj)
	{
		if (old_obj.length > 0)
		{
			length = old_obj.length;
			e = new int[length];
			for (unsigned int i = 0; i != length; i++)
			{
				e[i] = old_obj.e[i];
			}
		}
		else
		{
			e = nullptr;
			length = 0;
		}
	}

	~basisVector() 
	{
		if (length > 0)
			delete[] e;
	}

	/**
	* Compare operator
	* @param b BasisVector to compare to
	*/
	bool operator== (const basisVector &b) const
	{
		if (b.length != length)
			return false;
		else
		{
			for (unsigned int i = 0; i != length; i++)
			{
				if (e[i] != b.e[i])
					return false;
			}
		}
		return true;
	}

	/**
	* Assign operator
	* @param rhs BasisVector from which a copy is made
	*/
	basisVector& operator= (const basisVector &rhs)
	{

		if (rhs.length == length)
		{
			for (unsigned int i = 0; i != length; i++)
			{
				e[i] = rhs.e[i];
			}
		}
		else
		{
			if (e != nullptr)
				delete[] e;
			e = new int[rhs.length];
			length = rhs.length;
			for (unsigned int i = 0; i != length; i++)
			{
				e[i] = rhs.e[i];
			}
		}

		return *this;
	}


};


/**
* Hash class for vectors
*/
struct basisVectorHasher
{
	/**
	* Hash operator 
	* @param b Vector which is going to get hashed
	*/
	std::size_t operator()(const basisVector &b) const
	{
		std::size_t seed = 0;
		seed = boost::hash_range(b.e, b.e + b.length);
		return seed;
	}
};

/**
 * Represents a basis states, which consists of a basis vector and a complex coefficient.
 */
struct basisState
{
public:
	basisVector b;
	std::complex<double> coef;

	/**
 	* Creates a basis state by copying an existing one
 	* @param rhs The existing basis state to be copied
	 */
	basisState(const basisState &rhs)
	{
		b = rhs.b;
		coef = rhs.coef;
	}

	/**
 	* Creates a basis state of length 0
	 */
	basisState()
	{
		b = basisVector();
		coef = 0;
	}

	/**
 	* Creates a basis state by copying an existing one but changing its coefficient
 	* @param rhs The existing basis state to be copied
 	* @param c The new coefficient
	 */
	basisState(basisVector &rhs, std::complex<double> &c)
	{
		b = basisVector(rhs);
		coef = c;
	}

};

/**
 * Mother class for set of basis states
 */
class basicBasis
{
public:
    int numberModes;
    int qudits;
    int numberElements;
    basisVector* basisElements;
    std::unordered_map<basisVector, int, basisVectorHasher>  hashTable;
    
    
    virtual void createHashTable() = 0;
    basicBasis(){};
    virtual ~basicBasis(){};
    
};

/**
 * Contains a set of basis states
 */
class basis : public basicBasis
{
public:

	basis(int N, int K, int nbQudits, int capacity, int quant);
	basis(int N, int K, int nbQudits, int capacity);


	/**
 * Computes binomial coefficient (n; k)
 * @param n The first argument of the binomial coefficient
 * @param k The second argument of the binomial coefficient
 */
	int nchoosekSmart(int n, int k)
	{
		if (k > n)
			return 1;
		else
			return (int) boost::math::binomial_coefficient<double>(n, k);
	}


	basis& operator= (const basis& rhs)
	{
		if (numberElements > 0)
			delete[] basisElements;
		
		numberElements = rhs.numberElements;
		numberModes = rhs.numberModes;
		qudits = rhs.qudits;

		basisElements = new basisVector[numberElements];
		for (int i = 0; i != numberElements; i++)
			basisElements[i] = rhs.basisElements[i];
		return *this;
	}

	~basis();

     void createHashTable();

};

/**
 * Contains a set of basis states formed from the tensor product of two sets of basis states
 */
class tensorBasis : public basicBasis
{
public:

	tensorBasis(int N1, int K1, int N2, int K2, int C2);
    tensorBasis(){numberElements =0; numberModes = 0; qudits = -1; basisElements = nullptr;};
    
    void createHashTable();

	~tensorBasis()
	{
		if (numberElements != 0)
			delete[] basisElements;
	}

};

