#include "ellipticalharmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

EllipticalHarmonicOscillator::EllipticalHarmonicOscillator(System* system, double omega, double gamma) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_gamma   = gamma;
}

double EllipticalHarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* This function calculated the local energy for the wavefunction with particles 
    in the positions given by the input. The potential energy is always calculated the 
    same way, but the kinetic energy can be calcualted both analytically and numerically.
    This is determined by a bool statement parameter in the System class. */

    double hbar = 1.0;        // Planck's constant, but in natural units
    double m = 1.0;           // The boson mass, but in natural units

    double rSum2 = 0.0;
    double doubleDerivative = 0.0;

    double kineticEnergy;
    double potentialEnergy;

    particles = m_system->getParticles();

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <double> particlePosition(numberOfDimensions);

    for (int p1 = 0; p1 < numberOfParticles; p1++){
        particlePosition = particles[p1]->getPosition();

        for (int p2 = 0; p2 < numberOfDimensions; p2++){
            if (p2 == 2){ // This is for the elliptical trap case, else beta = 1
                rSum2 += m_gamma*particlePosition[p2]*particlePosition[p2];
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

    return potentialEnergy + kineticEnergy;
}

std::vector<double> EllipticalHarmonicOscillator::computeQuantumForce(std::vector<class Particle*> particles){
    /* This function calculates the quantum force/drift force with is used for importance
        sampling. The quantum force is given by the derivative of the wavefunction. */
    
     auto derivative = m_system->getWaveFunction()->computeDerivative(particles);

    for (int m=0;m<m_system->getNumberOfDimensions();m++){
        derivative[m] *= 2;
    }

    return derivative;
}

double EllipticalHarmonicOscillator::computeEnergyDerivative(std::vector<class Particle*> particles){
    /* OBS! This only calculates part of the derivative. The rest is done in Sampler */

    double expectationDerivative;

    expectationDerivative = computeLocalEnergy(particles)*m_system->getWaveFunction()->computeAlphaDerivative(particles);

    return expectationDerivative;
}