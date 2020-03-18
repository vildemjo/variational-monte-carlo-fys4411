#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double beta);

    double              computeLocalEnergy(std::vector<Particle*> particles);
    std::vector<double> computeQuantumForce(std::vector<class Particle*> particles);
    double              computeEnergyDerivative(std::vector<class Particle*> particles);
    double              getBeta(){ return m_beta; }

private:
    double m_omega = 0;
    double m_beta = 0;
};

