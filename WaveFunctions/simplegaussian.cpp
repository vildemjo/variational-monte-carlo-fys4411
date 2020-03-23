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

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {

    double rSum = 0.0;
    double beta = m_system->getHamiltonian()->getBeta();
    
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    particles = m_system->getParticles();

    for(int i1=0; i1<numberOfParticles; i1++){
        auto r = particles[i1]->getPosition();
        
        for(int n1=0; n1<numberOfDimensions; n1++){
            if (n1 == 2){
                rSum += beta*r[n1]*r[n1];
            }else{
                rSum += r[n1]*r[n1];
            }
        }
    }

    return exp(-m_parameters[0]*rSum);

}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
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
    double      rSum2               = 0.0;
    double      beta                = m_system->getHamiltonian()->getBeta();

    particles = m_system->getParticles();

    for(int i2=0; i2<numberOfParticles; i2++){
        auto r = particles[i2]->getPosition();
        for(int n2=0; n2<numberOfDimensions; n2++){
            if (n2 == 2){
                rSum2 += beta*r[n2]*r[n2];    
            }else{
                rSum2 += r[n2]*r[n2];
            }
        }   
    }
    
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2);
}

std::vector<double> SimpleGaussian::computeDerivative(std::vector<class Particle*> particles){
    
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector<double> vectorSum(numberOfDimensions);

    for (int i8 = 0; i8<m_system->getNumberOfParticles();i8++){

        auto r = particles[i8]->getPosition();
        double phiPart;
        double beta = m_system->getHamiltonian()->getBeta();

        for (int n8=0; n8<numberOfDimensions; n8++){
            if (n8 == 2){
                phiPart = -2*getParameters()[0]*beta*r[n8];
            }else{
                phiPart = -2*getParameters()[0]*r[n8];
            }
        vectorSum[n8] += phiPart;
        }
    }
    return vectorSum;


    

    
}

double SimpleGaussian::computeAlphaDerivative(std::vector<class Particle*> particles){

    double vectorSumSquared;
    double beta = m_system->getHamiltonian()->getBeta();

    for (int i10 = 0; i10<m_system->getNumberOfParticles();i10++){
        auto r = particles[i10]->getPosition();
        for (int n10=0; n10<m_system->getNumberOfDimensions(); n10++){
            if (n10 == 2){
                vectorSumSquared += beta*r[n10]*r[n10];
            }else{
                vectorSumSquared += r[n10]*r[n10];
            }
        }
    }

    return (-1)*vectorSumSquared;

}
