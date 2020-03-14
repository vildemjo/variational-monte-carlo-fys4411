#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {

    double hbar = 1.0;
    double m = 1.0;

    double rSum2 = 0.0;

    std::vector<double> r = std::vector<double>();
    
    for(int i2=0; i2<m_system->getNumberOfParticles(); i2++){
        r = particles[i2]->getPosition();
        for(int n2=0; n2<m_system->getNumberOfDimensions(); n2++){
            rSum2 += r[n2]*r[n2];
        }
    }

    double doubleDerivative;

    // std::cout << "evaluating rSum is ok" << std::endl;

    double potentialEnergy = 0.5*m*m_omega*rSum2;

    // Calculating the normalized second derivative either analytically or numerically
    if (m_system->getAnalytical() == true){
        doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
        // std::cout << "the double derivative is okay" << std::endl;
        }
    else{
        doubleDerivative = computeDoubleDerivativeNumerically(particles);
    }
    double kineticEnergy   = (-hbar*hbar/(2*m))*doubleDerivative;
    

    // std::cout << "no interaction energy calculation works" << std::endl;
    double interactionEnergy = m_system->getHamiltonian()->getInteractionPotential();

    return kineticEnergy + potentialEnergy + interactionEnergy;
}

std::vector<double> HarmonicOscillator::computeQuantumForce(std::vector<class Particle*> particles){

    /* Should maybe add the possibility of choosing nummerical derivative? */
    // std::cout << "into quantum F ok" << std::endl;
    /* computeDerivative checks if there is interaction or not and returns the relevant derivative */
    auto derivative = m_system->getWaveFunction()->computeDerivative(particles);

    // std::cout << "computeDerivative works" << std::endl;
    /* The quantum force is given by 2*psi_T^-1 derivative(psi_T) */ 
    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
    }

    return derivative;
}

double HarmonicOscillator::computeEnergyDerivative(std::vector<class Particle*> particles){
    /* OBS! This only calculates part of the derivative. The rest is done in Sampler */

    double expectationDerivative;

    expectationDerivative = computeLocalEnergy(particles)*m_system->getWaveFunction()->computeAlphaDerivative(particles);

    return expectationDerivative;
}
