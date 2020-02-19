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

    for(int i2=0; i2<numberOfParticles; i2++){
        for(int n2=0; n2<numberOfDimensions; n2++){
            rSum2 += r[i2][n2]*r[i2][n2];
        }
        
    }
    
    double potentialEnergy = 0.5*m*omega*rSum2;
    double kineticEnergy   = (-hbar*hbar/(2*m))*(-2*alpha*numberOfParticles*numberOfDimensions + 4*alpha*alpha*rSum2);
    return kineticEnergy + potentialEnergy;
}

