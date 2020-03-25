#pragma once
#include <vector>
#include <string>

class System {
public:
    bool metropolisStep             ();
    bool metropolisStepImportance   ();
    void runMetropolisSteps         (int numberOfMetropolisSteps, int firstCriteria, bool importanceOrNot, double stepLength);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setAnalytical              (bool statement);
    void setImportance              (bool statement);
    void setFileName                (std::string filename) {m_filename = filename;}
    void setHardCoreDiameter        (double hardCoreDiameter);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class InitialState*             getInitialState()   { return m_initialState; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    bool getAnalytical()                { return m_analytical; }
    bool getImportance()                { return m_importance; }
    double getStepLength()              { return m_stepLength; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getHardCoreDiameter()        { return m_hardCoreDiameter; }
    std::string getFileName()                { return m_filename; }
    double greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld);

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    bool                            m_analytical = false;
    bool                            m_importance = false;
    double                          m_diffConstant = 0.5;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_timeStep   = 0.1;
    double                          m_hardCoreDiameter = 0.5;
    std::string                     m_filename = "Output/no_name_specified";
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
};

