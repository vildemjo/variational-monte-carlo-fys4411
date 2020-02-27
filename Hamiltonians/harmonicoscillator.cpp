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
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double rSum2 = 0.0;

    std::vector<double> r = std::vector<double>();
    
    for(int i2=0; i2<m_system->getNumberOfParticles(); i2++){
        r = particles[i2]->getPosition();
        for(int n2=0; n2<m_system->getNumberOfDimensions(); n2++){
            rSum2 += r[n2]*r[n2];
        }
    }

    double doubleDerivative;

    double potentialEnergy = 0.5*m*omega*rSum2;

    // Calculating the normalized second derivative either analytically or numerically
    if (m_system->getAnalytical() == true){
        doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
        }
    else{
        doubleDerivative = computeDoubleDerivativeNumerically(particles);
    }
    double kineticEnergy   = (-hbar*hbar/(2*m))*doubleDerivative;
    return kineticEnergy + potentialEnergy;
}

std::vector<double> HarmonicOscillator::computeQuantumForce(std::vector<class Particle*> particles){

    /* Should maybe add the possibility of choosing nummerical derivative? */

    auto derivative = m_system->getWaveFunction()->computeDerivative(particles);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= -4*m_system->getWaveFunction()->getParameters()[0];
    }

    return derivative;
}

