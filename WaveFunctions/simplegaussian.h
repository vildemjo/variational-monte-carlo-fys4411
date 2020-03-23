#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);
    double computeAlphaDerivative(std::vector<class Particle*> particles);
    bool getDistanceCheck(){ return true; };
};
