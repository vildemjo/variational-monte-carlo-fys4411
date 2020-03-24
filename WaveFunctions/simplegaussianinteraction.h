#pragma once
#include "wavefunction.h"

class SimpleGaussianInteraction : public WaveFunction {
public:
    SimpleGaussianInteraction(class System* system, double alpha, double HardCoreDiameter);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);
    double computeAlphaDerivative(std::vector<class Particle*> particles);
    bool getDistanceCheck(std::vector<class Particle*> particles) { return calculateInterparticleDistances(particles); };

private:
    double evaluateCorrelationPart(std::vector<class Particle*> particles);
    double computeInteractionPartOfDoubleDerivative(std::vector<class Particle*> particles);
    std::vector <double> computeDerivativeOfu(std::vector<class Particle*> particles, int particleNumber);
    std::vector <double> computeDerivativeOneParticle(std::vector<class Particle*> particles, int particleIndex);
    bool calculateInterparticleDistances(std::vector<class Particle*> particles);
    void setDistances(std::vector<std::vector<double>> distances);
    void updateDistances(int particleNumber);
    std::vector<double> evaluateDifferenceVector();
    std::vector<std::vector<double>> getDistances(std::vector<class Particle*> particles) {
        calculateInterparticleDistances(particles);
        return m_distances; }
    std::vector<std::vector<double>> m_distances;
};
