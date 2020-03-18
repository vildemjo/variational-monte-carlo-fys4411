#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double beta) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_beta   = beta;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {

    double hbar = 1.0;
    double m = 1.0;

    double rSum2 = 0.0;
    double doubleDerivative = 0.0;

    double kineticEnergy;
    double potentialEnergy;

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double beta = m_system->getHamiltonian()->getBeta();

    std::vector <double> particlePosition(numberOfDimensions);

    for (int p1 = 0; p1 < numberOfParticles; p1++){
        particlePosition = particles[p1]->getPosition();

        for (int p2 = 0; p2 < numberOfDimensions; p2++){
            if (p2 == 2){
                rSum2 += beta*particlePosition[p2]*particlePosition[p2];
            }else{
                rSum2 += particlePosition[p2]*particlePosition[p2];
            }
        }
    }

    potentialEnergy = 0.5*m*m_omega*m_omega*rSum2;

    if (m_system->getAnalytical() == true){
        doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
    }else{
        doubleDerivative = computeDoubleDerivativeNumerically(particles);
    }

    kineticEnergy = -0.5*(hbar*hbar/m)*doubleDerivative;

    // std::vector<double> r = std::vector<double>();
    
    // for(int i2=0; i2<m_system->getNumberOfParticles(); i2++){
    //     r = particles[i2]->getPosition();
    //     for(int n2=0; n2<m_system->getNumberOfDimensions(); n2++){
    //         if (n2 == 2){
    //             rSum2 += m_beta*r[n2]*r[n2];
    //         }else{
    //             rSum2 += r[n2]*r[n2];
    //         }
    //     }
    // }

    // double doubleDerivative;

    // // std::cout << "evaluating rSum is ok" << std::endl;

    // double potentialEnergy = 0.5*m*m_omega*m_omega*rSum2;

    // // Calculating the normalized second derivative either analytically or numerically
    // if (m_system->getAnalytical() == true){
    //     doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
    //     // std::cout << "the double derivative is okay" << std::endl;
    //     }
    // else{
    //     doubleDerivative = computeDoubleDerivativeNumerically(particles);
    // }
    // double kineticEnergy   = (-hbar*hbar/(2.0*m))*doubleDerivative;
    
    // double interactionEnergy = 0;

    // // This should now not be neccecary because the steps where the Interaction potential
    // // is > 0, should not be accepted
    // if (m_system->getInteractionOrNot() == true){
    //     interactionEnergy = m_system->getHamiltonian()->getInteractionPotential();

    //     if (interactionEnergy > 0){
    //         std::cout << "interaction energy above zero: "<< interactionEnergy << std::endl;
    //     }
    // }

    return kineticEnergy + potentialEnergy; //+ interactionEnergy;
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
