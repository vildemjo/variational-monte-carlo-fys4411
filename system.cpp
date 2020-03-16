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

using namespace std;

bool System::metropolisStep() {
    
    // Evaluating the wave function for the present positions
    double oldWaveFunction = m_waveFunction->evaluate(m_particles);
    
    /* Here the type of sampling is checked 
    (Options: Brute Force Monte Carlo or Importance Sampling) */
    
    if (m_importanceSampling == false){  // Brute Force MC sampling
        
        std::vector<double> randomAmount = std::vector<double>();  // The distance the particle is moved in each direction

        // Choose a random particle
        int randomParticleIndex = Random::nextInt(m_numberOfParticles-1);


        // Change particle's position in all dimentions
        for(int m1=0;m1<m_numberOfDimensions; m1++){
            randomAmount.push_back(m_stepLength*(Random::nextDouble()-0.5));
            m_particles[randomParticleIndex]->adjustPosition(randomAmount[m1], m1);
        }

        // Calculate the wave function when the random particles is in it's new position
        double newWaveFunction = m_waveFunction->evaluate(m_particles);

        // Compare new wavefunction with old wavefunction
        if (Random::nextDouble() <= newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction)){
            // The step was accepted and will be counted in Sampler
            return true;
            }

        // Move the particle back to it's original position if not accepted
        for(int m2=0;m2<m_numberOfDimensions; m2++){
            m_particles[randomParticleIndex]->adjustPosition(-randomAmount[m2], m2);
        }
        // std::cout << "the step was not accepted" << std::endl;
        // The step was not accepted and not be counted
        return false;
    }
    else{ // Importance Sampling

        /* Here I extract the quantumForce for the old position r = (x,y,z) 
        and the new position r2 = (x2, y2, z2). I move the particle to the new position
        and then calculate both the wavefunction and the quantum force */

        std::vector<double> importanceAmount(m_numberOfDimensions); // The distance the particle is moved in each direction

        auto oldQuantumForce = m_hamiltonian->computeQuantumForce(m_particles);
        // std::cout << "quantum force works" << std::endl;

        // Choose a random particle
        int particleIndex = Random::nextInt(m_numberOfParticles-1);

        // Get the position of the particle
        auto oldPosition = m_particles[particleIndex]->getPosition();

        
        // Change position according to importance sampling:
        for(int m1=0;m1<m_numberOfDimensions; m1++){
            importanceAmount[m1] = (m_diffConstant*m_timeStep*oldQuantumForce[m1] + Random::nextGaussian(0, 1)*sqrt(m_timeStep));
            // std::cout << "importance amount works" << std::endl;
            m_particles[particleIndex]->adjustPosition(importanceAmount[m1], m1);
        }

        // Evaluate new wavefunction, quantum force and new position
        double newWaveFunction = m_waveFunction->evaluate(m_particles);
        std::vector<double>  newQuantumForce = m_hamiltonian->computeQuantumForce(m_particles);
        auto newPosition = m_particles[particleIndex]->getPosition();

        // Calculate the fraction of the Green's Functions for new to old and old to new
        double greensFunctionFrac = greensFunctionFraction(oldPosition, oldQuantumForce, newPosition, newQuantumForce);

        // Compare new wavefunction with old wavefunction including the Green's function fraction
        if (Random::nextDouble() <= greensFunctionFrac*newWaveFunction*newWaveFunction/(oldWaveFunction*oldWaveFunction)){
            // The move was accepted
            return true;
            }

        // Move back if it is not accepted
        for(int m4=0;m4<m_numberOfDimensions; m4++){
            m_particles[particleIndex]->adjustPosition(-importanceAmount[m4], m4);
        }
    // The move was not accepted
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

    // std::cout << "finished setting up" << std::endl;

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        
        // Counting the number of steps that are accepted by checking if it was accepted

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
    
    std::cout << "Out of MC loop" << std::endl;

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

void System::setInteractionOrNot(bool statement, double hardCoreDiameter){
    m_interactionOrNot = statement;
    m_hardCoreDiameter = hardCoreDiameter;
}

double System::greensFunctionFraction(std::vector<double> posNew, std::vector<double> posOld, std::vector<double> forceNew, std::vector<double> forceOld){
    double exponent;
    for (int n10=0; n10<m_numberOfDimensions; n10++){
        exponent += 0.5*(posNew[n10]*(forceNew[n10]-forceOld[n10]) + posOld[n10]*(forceNew[n10] + forceOld[n10]) + 0.5*m_diffConstant*m_timeStep*(forceOld[n10]*forceOld[n10]-forceNew[n10]*forceNew[n10]));
    }
    return exp(exponent);
}

