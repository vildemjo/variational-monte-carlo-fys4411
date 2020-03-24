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
#include <iostream>
#include <omp.h>

using namespace std;

bool System::metropolisStep() {
    /* This function contains the actual metropolis step. Here a random particle is
     chosen and moved randomly in the dimensions included in the simulation. Afterwards
     the acceptance of the move is checked using the standard Metropolis algorithm check. 
     If the move is accepted, the function returns true and lets function runMetropolisstep
     know that the step was accepted. If it is not accepted the particle is moved back to 
     its original position and the function return false. The implies that the sampling is 
     done with the wavefunction made up of particles at the original position.*/


    // A vector to save the distances the particle is moved in
    // case it has to be moved back
    std::vector<double> randomAmount = std::vector<double>();  

    double   oldWaveFunction      = m_waveFunction->evaluate(m_particles);
    int      randomParticleIndex  = Random::nextInt(m_numberOfParticles-1);

    
    for(int m1=0;m1<m_numberOfDimensions; m1++){
        randomAmount.push_back(m_stepLength*(Random::nextDouble()-0.5));
        m_particles[randomParticleIndex]->adjustPosition(randomAmount[m1], m1);
    }

    double newWaveFunction = m_waveFunction->evaluate(m_particles);
    
    
    if (Random::nextDouble() <= newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction)){
        return true;
        }

    for(int m2=0;m2<m_numberOfDimensions; m2++){
        m_particles[randomParticleIndex]->adjustPosition(-randomAmount[m2], m2);
    }
    
    return false;
}

bool System::metropolisStepImportance() {
    /* This function contains the actual metropolis step. Here a random particle is
     chosen and moved randomly in the dimensions included in the simulation. Afterwards
     the acceptance of the move is checked using importance sampling/Metropolis-Hastings
     algorithm check which includes Green's function and a quantum force/drift force 
     If the move is accepted, the function returns true and lets function 
     runMetropolisstepImportance know that the step was accepted. If it is not accepted the
     particle is moved back to its original position and the function return false. The 
     implies that the sampling is done with the wavefunction made up of particles at the 
     original position.*/


    // A vector to save the distances the particle is moved in
    // case it has to be moved back
    std::vector<double> importanceAmount(m_numberOfDimensions);

    int particleIndex = Random::nextInt(m_numberOfParticles-1);

    double oldWaveFunction    = m_waveFunction->evaluate(m_particles);
    auto   oldQuantumForce    = m_hamiltonian->computeQuantumForce(m_particles);
    auto   oldPosition        = m_particles[particleIndex]->getPosition();


    for(int m1=0;m1<m_numberOfDimensions; m1++){
        importanceAmount[m1] = (m_diffConstant*m_timeStep*oldQuantumForce[m1] 
                                + Random::nextGaussian(0, 1)*sqrt(m_timeStep));
        m_particles[particleIndex]->adjustPosition(importanceAmount[m1], m1);
    }

 
    double newWaveFunction  = m_waveFunction->evaluate(m_particles);
    auto   newQuantumForce  = m_hamiltonian->computeQuantumForce(m_particles);
    auto   newPosition      = m_particles[particleIndex]->getPosition();

    double greensFunctionFrac = greensFunctionFraction(oldPosition, oldQuantumForce,
                                                         newPosition, newQuantumForce);

    if (Random::nextDouble() <= greensFunctionFrac*newWaveFunction*newWaveFunction
                                            /(oldWaveFunction*oldWaveFunction)){
        return true;
        }

    
    for(int m4=0;m4<m_numberOfDimensions; m4++){
        m_particles[particleIndex]->adjustPosition(-importanceAmount[m4], m4);
    }

    return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps, int firstCriteria, double stepLength) {
    /* This function runs through the Monte Carlo cycles and performs the metropolis steps
    through the function metropolisStep. Here the energy and the information needed to evaluate
    the one-body density is sampled in the Sampler class and the result is printed to file. */
    
    m_particles                             = m_initialState->getParticles();
    m_sampler                               = new Sampler(this);
    m_numberOfMetropolisSteps               = numberOfMetropolisSteps;
    m_stepLength                            = stepLength;
    m_sampler->setNumberOfMetropolisSteps   (numberOfMetropolisSteps);
    m_sampler->setFileOutput                (firstCriteria);

    for (int i = 0; i < m_numberOfMetropolisSteps; i++) {
        
        bool acceptedStep = metropolisStep();

        // m_sampler->sample(acceptedStep);
        m_sampler->sampleAllEnergies(acceptedStep);
    }
    std::cout << "finished MC loop for alpha "<< getWaveFunction()->getParameters()[0] << std::endl;
    
    m_sampler->printOutputToEnergyFile();
    m_sampler->printOneBodyDensityToFile();
    m_sampler->computeAverages();
    // m_sampler->printOutputToEnergyAlphaFile();
}

void System::runMetropolisStepsImportance(int numberOfMetropolisSteps, int firstCriteria, double timeStep) {
    /* This function runs through the Monte Carlo cycles and performs the metropolis steps
    through the function metropolisStepImportance. Here the energy and the information needed 
    to evaluate the one-body density is sampled in the Sampler class and the result is printed 
    to file. */
    
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_timeStep                  = timeStep;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    m_sampler->setFileOutput(firstCriteria);

    
     for (int i=0; i < numberOfMetropolisSteps; i++) {
        
        bool acceptedStep = metropolisStepImportance();

        // m_sampler->sample(acceptedStep);
        m_sampler->sampleAllEnergies(acceptedStep);
    }
    
    m_sampler->printOutputToEnergyFile();
    m_sampler->printOneBodyDensityToFile();
    m_sampler->computeAverages();
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

double System::greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld){
    /* This function calculates the fraction between the Green's function for the transition
    from the old state to the new and the transition from the new to the old. This expression
    is used to determine wether a move is accepted or not when importance sampling is used. */
    
    double exponent;

    for (int n10=0; n10<m_numberOfDimensions; n10++){
        exponent += 0.5*(posNew[n10]*(forceNew[n10]-forceOld[n10]) 
                + posOld[n10]*(forceNew[n10] + forceOld[n10]) 
                + 0.5*m_diffConstant*m_timeStep
                *(forceOld[n10]*forceOld[n10]-forceNew[n10]*forceNew[n10]));
    }
    return exp(exponent);
}

void System::setHardCoreDiameter(double hardCoreDiameter){
    m_hardCoreDiameter = hardCoreDiameter;
}