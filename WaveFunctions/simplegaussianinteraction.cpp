#include "simplegaussianinteraction.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SimpleGaussianInteraction::SimpleGaussianInteraction(System* system, double alpha, double hardCoreDiameter) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_system->setHardCoreDiameter(hardCoreDiameter);
}

double SimpleGaussianInteraction::evaluate(std::vector<class Particle*> particles) {

    double rSum = 0.0;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();
    double beta = m_system->getHamiltonian()->getBeta();

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

    double interactionPart = evaluateCorrelationPart(particles);

    return exp(-m_parameters[0]*rSum)*interactionPart;

}

double SimpleGaussianInteraction::evaluateCorrelationPart(std::vector<class Particle*> particles) {

    int     numberOfParticles       = m_system->getNumberOfParticles();
    int     numberOfDimensions      = m_system->getNumberOfDimensions();
    int     uSumCheck               = 0;
    double  beta                    = m_system->getHamiltonian()->getBeta();
    double  correlationPart         = 1;
    double  uSum                    = 0;
    double  distances_j1_j2;
    auto    a                       = m_system->getHardCoreDiameter();
    auto    distances               = getDistances(particles);
    
    std::vector<double> distances_j1(numberOfParticles);
    std::vector<double> distances_j2(numberOfParticles);


    for (int j1 = 0; j1 < numberOfParticles-1; j1++){

        distances_j1 = distances[j1];

        for (int j2 = j1+1; j2 <numberOfParticles; j2++){
            distances_j1_j2 = distances_j1[j2];
            if ( distances_j1_j2 <= a ) {
                    uSumCheck += 1;
            }else{
                uSum += log(1-a/distances_j1_j2);
            }
        }
    }
    
    // Here the interaction part is set to zero directly (instead of exp(-infty))
    // if one of the distances are less than a
    if (uSumCheck > 0){
        correlationPart = 0;
    }else{
        correlationPart = exp(uSum);
    }
    // std::cout << "correlation part works" << std::endl;
    return correlationPart;

}

double SimpleGaussianInteraction::computeDoubleDerivative(std::vector<class Particle*> particles) {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double rSum2 = 0.0;
    double interactionPart = 0;
    double beta = m_system->getHamiltonian()->getBeta();

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

    interactionPart = computeInteractionPartOfDoubleDerivative(particles);
    
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2) + interactionPart;
}

std::vector<double> SimpleGaussianInteraction::computeDerivative(std::vector<class Particle*> particles){
    
    int                     numberOfDimensions          = m_system->getNumberOfDimensions();
    double                  phiPart                     = 0;
    double                  beta                        = m_system->getHamiltonian()->getBeta();
    std::vector <double>    vectorWithInteraction       (numberOfDimensions);
    
    
    for (int i8 = 0; i8<m_system->getNumberOfParticles();i8++){

        auto uDerivative = computeDerivativeOfu(particles, i8);
        auto r = particles[i8]->getPosition();

        for (int n8=0; n8<numberOfDimensions; n8++){
            if (n8 == 2){
                phiPart = -2*getParameters()[0]*r[n8]*beta;
            }else{
                phiPart = -2*getParameters()[0]*r[n8];
            }
            vectorWithInteraction[n8] += phiPart + uDerivative[n8];
        }
    }

    return vectorWithInteraction;
}

double SimpleGaussianInteraction::computeAlphaDerivative(std::vector<class Particle*> particles){

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

double SimpleGaussianInteraction::computeInteractionPartOfDoubleDerivative(std::vector<class Particle*> particles){

    double      firstTerm, secondTerm, thirdTerm;
    int         numberOfParticles = m_system->getNumberOfParticles();
    int         numberOfDimentions = m_system->getNumberOfDimensions();
    double      uDerivative = 1;
    double      uDoubleDerivative = 1;


    for(int l5 = 0; l5 < numberOfParticles; l5++){ 

        auto uPart = computeDerivativeOfu(particles, l5);
        auto derivativePhi = computeDerivativeOneParticle(particles, l5);

        for (int l4 = 0; l4<numberOfDimentions; l4++){
            firstTerm += derivativePhi[l4]*uPart[l4];
            secondTerm += uPart[l4]*uPart[l4];
        }
        thirdTerm += uPart[-1];
        
    }

    return 2*firstTerm + secondTerm + thirdTerm;
}

std::vector <double> SimpleGaussianInteraction::computeDerivativeOfu(std::vector<class Particle*> particles, int particleNumber){

    int                     numberOfParticles       = m_system->getNumberOfParticles();
    int                     numberOfDimentions      = m_system->getNumberOfDimensions();
    double                  uDerivative             = 0;
    double                  a                       = m_system->getHardCoreDiameter();
    double                  uDoubleDerivative       = 0;
    double                  uTotalDoubleDerivative  = 0;
    std::vector <double>    uTotalDerivative        (numberOfDimentions);
    std::vector <double>    uAllStuff               (numberOfDimentions+1);

    auto ri         = particles[particleNumber]->getPosition();                              // (x_i, y_i, z_i)
    auto difference = getDistances(particles)[particleNumber];
    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        if (particleNumber != l1){
            auto rj         = particles[l1]->getPosition();                          // (x_j, y_j, z_j)
            auto rLength    = difference[l1];                                   // r_ij
    
            /* Here sum u'(r_ij) is determined based on the relationship 
            between r_ij (distance between particles) and a (hard core diameter) */

            if (rLength <= a){
                uDerivative = -1e50;
                uDoubleDerivative = -1e50;
                std::cout << "r_ij < a in computeDerivativeOfu" << std::endl;
            }else{
                uDerivative = -a/(a*rLength-rLength*rLength);
                uDoubleDerivative = a*(a-2*rLength)/(rLength*rLength*(a-rLength)*(a-rLength));
            }

            for (int l3 = 0; l3<numberOfDimentions; l3++){
                uTotalDerivative[l3] += ((ri[l3]-rj[l3])/rLength)*uDerivative;
            }

            uTotalDoubleDerivative += uDoubleDerivative + ((numberOfDimentions-1)/rLength)*uDerivative;
        }
    }

    for (int l4 = 0; l4<numberOfDimentions+1; l4++){
        if(l4<numberOfDimentions){
            uAllStuff[l4] = uTotalDerivative[l4];
        }else{
            uAllStuff[l4] = uTotalDoubleDerivative;
        }
    }

    return uAllStuff;
}

std::vector<double> SimpleGaussianInteraction::computeDerivativeOneParticle(std::vector<class Particle*> particles, int particleIndex){
    
    int                 numberOfDimensions          = m_system->getNumberOfDimensions();
    double              beta                        = m_system->getHamiltonian()->getBeta();
    std::vector<double> derivativeVector            (numberOfDimensions);

    auto r = particles[particleIndex]->getPosition();

    for (int j3=0; j3<numberOfDimensions; j3++){
        if (j3 == 2){
            derivativeVector[j3] = -2*getParameters()[0]*beta*r[j3];
        }else{
            derivativeVector[j3] = -2*getParameters()[0]*r[j3];
        }
    }

    return derivativeVector;
    
}

bool SimpleGaussianInteraction::calculateInterparticleDistances(std::vector<class Particle*> particles){
    std::vector <std::vector <double>>                distances       (m_system->getNumberOfParticles()); // a matrix of the distance between all particles 
    std::vector <double>                              difference      (m_system->getNumberOfParticles()); // the distance between particle j and all other particles i where j>i
    std::vector <double>                              vectorDistance  (m_system->getNumberOfDimensions()); 

    double a = m_system->getHardCoreDiameter();
    m_system->getHamiltonian()->setInteractionPotential(false);
    int distanceCheck = 0;
    std::vector <double> r1(m_system->getNumberOfDimensions()), r2(m_system->getNumberOfDimensions());

    for (int k1 = 0; k1 < m_system->getNumberOfParticles(); k1++){
        // std::cout << "entering loop. Particle 1 is " << k1 << std::endl;

        r1 = particles[k1]->getPosition();
        // std::cout << "got position" << std::endl;
        for (int k2 = 0; k2 <m_system->getNumberOfParticles(); k2++){
            // std::cout << "entering loop. Particle 2 is " << k2 << std::endl;
            // std::cout << "hence calculating: r_" << k1 << k2 << std::endl;
            if (k1 != k2){
                r2 = particles[k2]->getPosition();
                for (int k3 = 0; k3<m_system->getNumberOfDimensions(); k3++){
                    difference[k2] +=  (r1[k3]-r2[k3])*(r1[k3]-r2[k3]);              // (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2
                }
                difference[k2] = sqrt(difference[k2]);                                  // sqrt((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2)

                // std::cout << "this is r_ij for p " << k1+1 << std::endl;
                // std::cout << difference[0] << ",\t" << difference[1] << ",\t" << std::endl;
        

                if (difference[k2] < a){
                    std::cout << "too small difference: " << difference[k2] << "<" << a << std::endl;
                    distanceCheck +=1;
                }
            }
        }
        // std::cout << "all distances are big enough" << std::endl;

        distances[k1] = difference;

        // std::cout << "have saved the vector of distance to the matrix of distances ones" << std::endl;

    }
    setDistances(distances);

    if (distanceCheck > 0){
        m_system->getHamiltonian()->setInteractionPotential(true);  // Telling the interaction potential that a distance is smaller than a
        return false;
    }
    // std::cout << "have set the distances and they are bigger than a" << std::endl;
    return true;
}


void SimpleGaussianInteraction::setDistances(std::vector<std::vector<double>> distances){
    m_distances = distances;
    // std::cout << "works inside setDistance function too" << std::endl;
}


std::vector<double> SimpleGaussianInteraction::evaluateDifferenceVector(){

    std::vector <double> differenceVector(m_system->getNumberOfDimensions());

    for (int k1 = 0; k1 < m_system->getNumberOfParticles()-1; k1++){
        auto r1 = m_system->getParticles()[k1]->getPosition();
        for (int k2 = k1+1; k2 <m_system->getNumberOfParticles(); k2++){
            auto r2 = m_system->getParticles()[k2]->getPosition();
            for (int k3 = 0; k3<m_system->getNumberOfDimensions(); k3++){
                differenceVector[k3] += r1[k3]-r2[k3];
            }
        }
    }
    return differenceVector;
}