#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void printOutputToFile();
    void setFileOutput(int firstCriteria);
    void computeAverages();
    double getEnergy()          { return m_energy; }
    double getDerivative()      { return m_derivative; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_energySquared = 0;
    double  m_derivative = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_cumulativeEnergyDerivative = 0;
    double  m_cumulativeAlphaDerivative = 0;
    int     m_numberOfAcceptedSteps = 0;
    int     m_firstCriteria = 1;
    // string  m_filename;
    bool    m_printOrNot;
    class System* m_system = nullptr;
};
