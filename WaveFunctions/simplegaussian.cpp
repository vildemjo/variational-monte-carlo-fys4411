#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particles[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    double rSum = 0.0;

    int numberOfParticles = particles.size();
    int numberOfDimensions = particles[0]->getPosition().size();

    // std::vector<double> r = std::vector<double>();

    for(int i1=0; i1<numberOfParticles; i1++){
        auto r = particles[i1]->getPosition();
        for(int n1=0; n1<numberOfDimensions; n1++){
            rSum += r[n1]*r[n1];
        }
    }
    return exp(-m_parameters[0]*rSum);

    return 0;
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

    int numberOfParticles = particles.size();
    int numberOfDimensions = particles[0]->getPosition().size();

    double rSum2 = 0.0;

    for(int i2=0; i2<numberOfParticles; i2++){
        auto r = particles[i2]->getPosition();
        for(int n2=0; n2<numberOfDimensions; n2++){
            rSum2 += r[n2]*r[n2];
        }   
    }
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2);
}

std::vector<double> SimpleGaussian::computeDerivative(std::vector<class Particle*> particles){
    
    std::vector<double> vectorSum(m_system->getNumberOfDimensions());

    for (int i8 = 0; i8<m_system->getNumberOfParticles();i8++){
        auto r = particles[i8]->getPosition();
        for (int n8=0; n8<m_system->getNumberOfDimensions(); n8++){
            vectorSum[n8] += r[n8];
        }
    }

    return vectorSum;
}


double SimpleGaussian::computeAlphaDerivative(std::vector<class Particle*> particles){

    double vectorSumSquared;

    for (int i10 = 0; i10<m_system->getNumberOfParticles();i10++){
        auto r = particles[i10]->getPosition();
        for (int n10=0; n10<m_system->getNumberOfDimensions(); n10++){
            vectorSumSquared += r[n10]*r[n10];
        }
    }

    return (-1)*vectorSumSquared;

}
