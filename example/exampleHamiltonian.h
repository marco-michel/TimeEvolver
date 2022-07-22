#pragma once

#include <vector>
#include <complex>

#include "hamiltonian.h"
#include "Basis.h"

/**
 * The specific exemplary Hamiltonian as defined in the accompnying paper
 */ 
class exampleHamiltonian : public Hamiltonian
{
public:
    int C; int N0; int Nm; double DeltaN; int Q1; int Q2; double C0; double CGap1; double CGap2; double Cm; double CmS;    
   
    exampleHamiltonian(int n0, int nm, double deltan, int q1, int q2, int c, double c0, double cgap1, double cgap2, double cm, double cms);
    void createSimplifiedHamiltonian();
  
   	basisVector createInitState();

private:
    double symBreak(int a, int b);
    std::vector<opTerm> simplifiedSphericalInteraction(int Q1, int Q2, double gm, double gms, int C);
};






