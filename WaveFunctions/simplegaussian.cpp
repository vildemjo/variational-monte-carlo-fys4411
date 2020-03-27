#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate() {

    auto rSum2 = calculatePositionSumSquared();

    return exp(-m_parameters[0]*rSum2);

}

double SimpleGaussian::computeDoubleDerivative() {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    int         numberOfParticles   = m_system->getNumberOfParticles();
    int         numberOfDimensions  = m_system->getNumberOfDimensions();

    auto rSum2 = calculatePositionSumSquared(); 
    
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2);
}

std::vector<double> SimpleGaussian::computeDerivative(int particleIndex){
    
    int numberOfDimensions = m_system->getNumberOfDimensions();
    auto m_particles = m_system->getParticles();

    std::vector<double> vectorSum(numberOfDimensions);

    auto r = m_particles[particleIndex]->getPosition();

    for (int n8=0; n8<numberOfDimensions; n8++){
        vectorSum[n8] = -2*getParameters()[0]*r[n8];
    }

    return vectorSum;

}

double SimpleGaussian::computeAlphaDerivative(){

    auto vectorSumSquared = calculatePositionSumSquared();

    return (-1)*vectorSumSquared;

}
