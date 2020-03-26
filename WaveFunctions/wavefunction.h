#pragma once
#include <vector>
#include "../system.h"


class WaveFunction {
public:
    WaveFunction(class System* system);
    void updateOneBodyDensity(std::vector<class Particle*> particles);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual std::vector<double> computeDerivative(int particleIndex, std::vector<class Particle*> particles) = 0;
    virtual double computeAlphaDerivative(std::vector<class Particle*> particles) = 0;
    virtual bool getDistanceCheck(std::vector<class Particle*> particles) = 0;
    std::vector <std::vector<double> > getOneBodyDensity(){ return m_oneBodyDensity; };
    void setOneBodyDensityBins(int numberOfBins, double densityLength);
    
protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;

private:
    std::vector<std::vector<double> > m_oneBodyDensity = std::vector<std::vector<double> >(); 
    int m_numberOfBins = 50;
    double m_densityLength = 1.0;
};

