#pragma once
#include <vector>


class System {
public:
    bool metropolisStep             ();
    void runMetropolisSteps         (int numberOfMetropolisSteps, bool statement, int firstCriteria);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    void setAnalytical              (bool statement);
    void setImportanceSampling      (bool statement, double timeStep);
    void setInteractionOrNot        (bool statement, double hardCoreDiameter);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    class InitialState*             getInitialState()   { return m_initialState; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    bool getAnalytical()                { return m_analytical; }
    bool getImportanceSampling()        { return m_importanceSampling; }
    bool getInteractionOrNot()          { return m_interactionOrNot; }
    double getStepLength()              { return m_stepLength; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getHardCoreDiameter()        { return m_hardCoreDiameter; }
    double greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld);

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    bool                            m_analytical = false;
    bool                            m_importanceSampling = false;
    bool                            m_printToFileOrNot = false;
    bool                            m_interactionOrNot = false;
    double                          m_diffConstant = 0.5;
    double                          m_timeStep = 0;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_hardCoreDiameter = 0.5;
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
};

