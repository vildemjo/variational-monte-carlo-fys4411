#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    
    // Calculate old wavefunction
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);

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

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        
        bool acceptedStep = metropolisStep();

        // Calculate energy (with choise of num or analytical?)

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        m_sampler->sample(acceptedStep);
    }
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
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


