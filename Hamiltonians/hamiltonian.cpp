#include "hamiltonian.h"
#include "../particle.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include <iostream>

using std::cout;
using std::endl;

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeDoubleDerivativeNumerically(std::vector<class Particle*> particles) {

    double step = 1e-3;
    double dpsidr2 = 0.0;
    double rSum2 = 0.0;

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double waveNext = 0;
    double waveLast = 0;
    double waveCurrent = 0;

    // Changing the copies to obtain an array of "next step" and "previous step"
    for(int i4=0; i4<numberOfParticles; i4++){
         for(int n4=0; n4<numberOfDimensions;n4++){
            // Wave function at forward step
            particles[i4]->adjustPosition(step, n4);
            waveNext = m_system->getWaveFunction()->evaluate(particles);

            // Wave function at backward step
            particles[i4]->adjustPosition(-2*step, n4);
            waveLast = m_system->getWaveFunction()->evaluate(particles);

            // Calculating the part of the double derivative which involves psi(x+dx) and psi(x-dx)
            dpsidr2 += (waveNext+waveLast)/(step*step);

            // Resetting the particlepositions so that a new particle and spesific dimension can be calculated
            particles[i4]->adjustPosition(step, n4); 
        }
    }

    waveCurrent = m_system->getWaveFunction()->evaluate(particles);

    // Calculating the part of the double derivative which involves psi(x)
    dpsidr2 += -2*numberOfParticles*numberOfDimensions*waveCurrent/(step*step);

    return (1/waveCurrent)*dpsidr2;
 
}

void Hamiltonian::setInteractionPotential(bool statement){
    if (statement == true){
        m_interactionPotential = 1e20;
    }else{
        m_interactionPotential = 0;
    }
}