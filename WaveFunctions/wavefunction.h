#pragma once
#include <vector>
#include "../system.h"


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeDerivative(std::vector<class Particle*> particles) = 0;
    virtual double computeAlphaDerivative(std::vector<class Particle*> particles) = 0;
    virtual bool getDistanceCheck() = 0;
    double computeInteractionPartOfDoubleDerivative(std::vector<class Particle*> particles);
    std::vector <double> computeDerivativeOfu(std::vector<class Particle*> particles, int particleNumber);
    std::vector <double> computeDerivativeOneParticle(std::vector<class Particle*> particles, int particleIndex);
    
protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

