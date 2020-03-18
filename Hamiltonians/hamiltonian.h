#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);

    double                      computeDoubleDerivativeNumerically(std::vector<class Particle*> particles);
    virtual double              computeLocalEnergy(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeQuantumForce(std::vector<class Particle*> particles) = 0;
    virtual double              computeEnergyDerivative(std::vector<class Particle*> particles) = 0;
    void                        setInteractionPotential(bool statement);
    double                      getInteractionPotential() { return m_interactionPotential; }
    virtual double              getBeta() = 0;

protected:
    class System* m_system = nullptr;
    int m_interactionPotential = 0;
};

