#include "exampleHamiltonian.h"
#include "Basis.h"

/**
 * Initializes the Hamiltonian for a specific choice of parameters (see accompanying paper for definition of parameters). This class allows for slightly more general choices of parameters than the Hamiltonian displayed in the paper.
 * @param n0 Parameter N0
 * @param nm Parameter Nm
 * @param deltan Parameter Delta N
 * @param q1 Parameter K
 * @param q2 Parameter K'
 * @param c Maximal occupation of quidt-modes; set it to 1
 * @param c0 Parameter C0
 * @param cgap1 Ratio of parameter epsilonm and square root of N0; usually it is 1
 * @param cgap2 Ratio of parameter epsilonm and square root of N0; usually it is 1
 * @param cm Parameter Cm
 * @param cms Parameter Cm
 */ 
exampleHamiltonian::exampleHamiltonian(int n0, int nm, double deltan, int q1, int q2, int c, double c0, double cgap1, double cgap2, double cm, double cms)
{
 
	C = c; N0 = n0; Nm = nm; DeltaN = deltan; Q1 = q1; Q2 = q2;
	C0 = c0; CGap1 = cgap1; CGap2 = cgap2; Cm = cm; CmS = cms; 
	Hamiltonian();
}

	/**
 * Builds the operators that make up the Hamiltonian
 * @return All operators which sum up to the Hamiltonian
 */ 
std::vector<opTerm> exampleHamiltonian::createSimplifiedHamiltonian()
{
    std::vector<opTerm> ret;
    opTerm tmp;
    
    std::vector<opTerm> ham;
    
    for(int i = 2; i < Q1+2; i++)
        ham.push_back(createNumberOperator(i, CGap1 * std::sqrt(N0)));
    
    for(int i = 2+Q1; i < Q1+Q2+2; i++)
        ham.push_back(createNumberOperator(i, CGap2 * std::sqrt(N0)));

    ham.push_back(linInteraction(0,1,0,0,true,C0));
    ham.push_back(linInteraction(0,1,0,0,false,C0));
    
    for(int i = 2; i < 2+Q1; i++)
        ham.push_back(createFourPointA(0,0,i,i,0,0,C,C,-std::sqrt(N0)*CGap1/N0));
    
    for(int i = 2+Q1; i < 2+Q1+Q2; i++)
        ham.push_back(createFourPointA(0,0,i,i,0,0,C,C,-std::sqrt(N0)*CGap2/(N0-DeltaN)));

    std::vector<opTerm> interactions = simplifiedSphericalInteraction(Q1, Q2, Cm, CmS, C);
    ham.insert(ham.end(),interactions.begin(),interactions.end());
    
    hamiltonOperator = ham;
    hamiltonianString = ham;
    return ham;
    
}

std::vector<opTerm> exampleHamiltonian::simplifiedSphericalInteraction(int Q1, int Q2, double gm, double gms, int C)
{
    std::vector<opTerm> terms;
    opTerm tmp;
    double coeff = 0;
    
    for(int m = 2; m <= Q1+Q2+1; m++)
    {
        for(int M = 2; M <= Q1+Q2+1; M++)
        {
            coeff = 0;
            if(m <= 1 + Q1 && M >Q1+1)
                coeff = symBreak(m,M)*gm;
            else if(((m<=Q1+1 && M <= 1+Q1) || (m>1+Q1 && M > 1+Q1)) && m<M)
                coeff = symBreak(m,M)*gms;
            if(coeff != 0)
            {
                tmp.coef.real(coeff);
                tmp.capacity.push_back(C);
                tmp.capacity.push_back(C);
                tmp.mode.push_back(m);
                tmp.mode.push_back(M);
                tmp.operations.push_back(-1);
                tmp.operations.push_back(1);
                
                terms.push_back(tmp);
                tmp.reset();
                
                tmp.coef.real(coeff);
                tmp.capacity.push_back(C);
                tmp.capacity.push_back(C);
                tmp.mode.push_back(M);
                tmp.mode.push_back(m);
                tmp.operations.push_back(-1);
                tmp.operations.push_back(1);
                
                terms.push_back(tmp);
                tmp.reset();
                
            }
        }
    }
    
    return terms;
}

double exampleHamiltonian::symBreak(int a, int b)
{
    double coef = std::fmod(std::sqrt(2)*std::pow(a,3)+std::sqrt(7)*std::pow(b,5),1.0);
    if (coef < 0.5)
        coef -= 1.0;
    return coef;
}

	/**
 * Creates a specific initial state, in which the first mode has occupation N0 and the first Nm qudit-modes have occupation 1
 * @return The initial state
 */ 
basisVector exampleHamiltonian::createInitState()
{
	
    basisVector ret;
    
   
	ret.length = Q1 + Q2 + 2;
	ret.e = new int[ret.length];

	int nbModes = Q1 + Q2 + 2;

	int minQNm = std::min(Q1, Nm);

	ret.e[0] = N0;
	ret.e[1] = 0;
	for (int i = 0; i != minQNm; i++)
		ret.e[i + 2] = 1;
	for (int i = 0; i != nbModes - 2 - minQNm; i++)
		ret.e[i + 2 + minQNm] = 0;


	return ret;
}
