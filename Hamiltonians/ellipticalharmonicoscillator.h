#pragma once
#include "hamiltonian.h"
#include <vector>

class EllipticalHarmonicOscillator : public Hamiltonian {
public:
    EllipticalHarmonicOscillator(System* system, double omega, double gamma);

    double              computeLocalEnergy(std::vector<Particle*> particles);
    std::vector<double> computeQuantumForce(std::vector<class Particle*> particles);
    double              computeEnergyDerivative(std::vector<class Particle*> particles);

private:
    double m_omega = 0;
    double m_gamma = 0;
};
