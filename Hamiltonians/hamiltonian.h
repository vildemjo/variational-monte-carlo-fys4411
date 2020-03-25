#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);

    double                      computeDoubleDerivativeNumerically(std::vector<class Particle*> particles);
    virtual double              computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeQuantumForce(std::vector<class Particle*> particles) = 0;



protected:
    class System* m_system = nullptr;
};

