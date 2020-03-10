#include "system.h"
#include <cassert>
#include <cmath>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <string>

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    
    // Calculate old wavefunction
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);
    
    if (m_importanceSampling == false){  // brute force MC sampling
        // Choose random particle
        int randomParticleIndex = Random::nextInt(m_numberOfParticles-1);

        std::vector<double> randomAmount = std::vector<double>();

        // Change particle's position in all dimentions
        for(int m1=0;m1<m_numberOfDimensions; m1++){
            randomAmount.push_back(m_stepLength*(Random::nextDouble()-0.5));
            m_particles[randomParticleIndex]->adjustPosition(randomAmount[m1], m1);
        }

        // Calculate new wavefunction
        double newWaveFunction = m_waveFunction->evaluate(m_particles);

        // Compare new wavefunction with old wavefunction
        if (Random::nextDouble() <= newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction)){
            return true;
            }

        // Move the particle back if not accepted
        for(int m2=0;m2<m_numberOfDimensions; m2++){
            m_particles[randomParticleIndex]->adjustPosition(-randomAmount[m2], m2);
        }

        return false;
    }
    else{ // Importance Sampling

        /* Here I have to extract the quantumForce for the old position r = (x,y,z) 
        and the new position r2 = (x2, y2, z2). I have to move the particle to the new position
        and then calculate both the wavefunction and the quantum force */

        std::vector<double> oldQuantumForce = m_hamiltonian->computeQuantumForce(m_particles);
        
        // Change position according to importance sampling stuff!
        int particleIndex = Random::nextInt(m_numberOfParticles-1);

        auto oldPosition = m_particles[particleIndex]->getPosition();

        std::vector<double> importanceAmount = std::vector<double>();

        // Change particle's position in all dimentions
        for(int m1=0;m1<m_numberOfDimensions; m1++){
            importanceAmount.push_back(m_diffConstant*m_timeStep*oldQuantumForce[m1] + Random::nextGaussian(0, 1)*sqrt(m_timeStep));
            m_particles[particleIndex]->adjustPosition(importanceAmount[m1], m1);
        }

        // Calculate new wavefunction and quantum force
        double newWaveFunction = m_waveFunction->evaluate(m_particles);
        std::vector<double>  newQuantumForce = m_hamiltonian->computeQuantumForce(m_particles);
        auto newPosition = m_particles[particleIndex]->getPosition();

        double greensFunctionFrac = greensFunctionFraction(oldPosition, oldQuantumForce, newPosition, newQuantumForce);

        // Compare new wavefunction with old wavefunction
        if (Random::nextDouble() <= greensFunctionFrac*newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction)){
            return true;
            }
        // Move back if it is not accepted
        for(int m4=0;m4<m_numberOfDimensions; m4++){
            m_particles[particleIndex]->adjustPosition(-importanceAmount[m4], m4);
        }
    
    return false;

    }
    
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, bool statement, int firstCriteria) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_sampler->setFileOutput(firstCriteria);
    m_printToFileOrNot = statement;

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
        * the m_sampler instance of the Sampler class. Make sure, though,
        * to only begin sampling after you have let the system equilibrate
        * for a while. You may handle this using the fraction of steps which
        * are equilibration steps; m_equilibrationFraction.
        */
        
        m_sampler->sample(acceptedStep);
    }
    m_sampler->computeAverages();
    
    if (m_printToFileOrNot == false){
    m_sampler->printOutputToTerminal();
    }
    else{
    m_sampler->printOutputToFile();
    }
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setAnalytical(bool statement){
    m_analytical = statement;
}

void System::setImportanceSampling(bool statement, double timeStep){
    m_importanceSampling = statement;
    m_timeStep = timeStep;
}

double System::greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld){
    double exponent;
    for (int n10=0; n10<m_numberOfDimensions; n10++){
        exponent += 0.5*(posNew[n10]*(forceNew[n10]-forceOld[n10]) + posOld[n10]*(forceNew[n10] + forceOld[n10]) + 0.5*m_diffConstant*m_timeStep*(forceOld[n10]*forceOld[n10]-forceNew[n10]*forceNew[n10]));
    }
    return exp(exponent);
}

