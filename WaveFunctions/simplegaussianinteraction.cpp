#include "simplegaussianinteraction.h"
#include "ellipticalharmonicoscillator.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "../particle.h"
#include "InitialStates/initialstate.h"
#include <iostream>

SimpleGaussianInteraction::SimpleGaussianInteraction(System* system, double alpha, double hardCoreDiameter, double beta) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
    m_system->setHardCoreDiameter(hardCoreDiameter);
    m_beta = beta;
}

double SimpleGaussianInteraction::evaluate() {

    auto rSum = calculatePositionSumSquared();

    double interactionPart = evaluateCorrelationPart();

    if (interactionPart == 0){
        // std::cout << "interactionpart is 0 in evaluate" << std::endl;
    }

    return exp(-m_parameters[0]*rSum)*interactionPart;

}

double SimpleGaussianInteraction::evaluateCorrelationPart() {

    int     numberOfParticles       = m_system->getNumberOfParticles();
    int     uSumCheck               = 0;
    double  correlationPart         = 1;
    double  uSum                    = 0;
    double  distances_j1_j2         = 0;
    auto    a                       = m_system->getHardCoreDiameter();
    auto    distances               = getDistances(m_system->getParticles());
    
    std::vector<double> distances_j1(numberOfParticles);
    std::vector<double> distances_j2(numberOfParticles);


    for (int j1 = 0; j1 < numberOfParticles-1; j1++){

        distances_j1 = distances[j1];

        for (int j2 = j1+1; j2 <numberOfParticles; j2++){
            distances_j1_j2 = distances_j1[j2];
            std::cout << "distance: " << distances_j1_j2 << std::endl;
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
        // std::cout << "interaction part set to 0" << std::endl;
    }else{
        correlationPart = exp(uSum);
    }
    // std::cout << "correlation part works" << std::endl;
    return correlationPart;

}

double SimpleGaussianInteraction::computeDoubleDerivative() {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double interactionPart = 0;

    auto rSum2 = calculatePositionSumSquared();

    interactionPart = computeInteractionPartOfDoubleDerivative();
    
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2) + interactionPart;
}

std::vector<double> SimpleGaussianInteraction::computeDerivative(int particleIndex){
    
    int                     numberOfDimensions          = m_system->getNumberOfDimensions();
    double                  phiPart                     = 0;
    std::vector <double>    vectorWithInteraction       (numberOfDimensions);
    std::vector <double>    uTotalDerivative            (numberOfDimensions);
    double uDerivative                                  = 0;
    double a = m_system->getHardCoreDiameter();

    auto m_particles = m_system->getParticles();

    auto ri = m_particles[particleIndex]->getPosition();

    auto theDifferences = getDistances(m_particles)[particleIndex];

    for (int n9 = 0; n9< (int) theDifferences.size(); n9++){
        auto rLength = theDifferences[n9];
        auto rj = m_particles[n9]->getPosition();

        if (n9 != particleIndex){
            if (rLength <= a){
                uDerivative = -1e50;
                std::cout << "r_ij < a in computeDerivative" << std::endl;
            }else{
                uDerivative = -a/(a*rLength-rLength*rLength);
            }

            for (int l3 = 0; l3<numberOfDimensions; l3++){
                uTotalDerivative[l3] += ((ri[l3]-rj[l3])/rLength)*uDerivative;
            }
        }
    }

    for (int n8=0; n8<numberOfDimensions; n8++){
        if (n8 == 2){
            phiPart = -2*getParameters()[0]*ri[n8]*m_beta;
        }else{
            phiPart = -2*getParameters()[0]*ri[n8];
        }
        vectorWithInteraction[n8] = phiPart + uTotalDerivative[n8];
    }

    return vectorWithInteraction;
}

double SimpleGaussianInteraction::computeAlphaDerivative(){
    auto m_particles = m_system->getParticles();

    auto vectorSumSquared =  calculatePositionSumSquared();

    return (-1)*vectorSumSquared; // No interaction part because it is divided away.

}


double SimpleGaussianInteraction::computeInteractionPartOfDoubleDerivative(){

    double      firstTerm = 0; double secondTerm = 0; double thirdTerm = 0;
    int         numberOfParticles = m_system->getNumberOfParticles();
    int         numberOfDimentions = m_system->getNumberOfDimensions();
    auto        m_particles = m_system->getParticles();


    for(int l5 = 0; l5 < numberOfParticles; l5++){ 

        auto derivativePhi = computeDerivativeOneParticle(l5);
        auto uPart = computeDerivativeOfu(l5);

        for (int l4 = 0; l4<numberOfDimentions; l4++){
            firstTerm += derivativePhi[l4]*uPart[l4];
            secondTerm += uPart[l4]*uPart[l4];
        }
        thirdTerm += uPart[-1];
        
    }

    return 2*firstTerm + secondTerm + thirdTerm;
}

std::vector <double> SimpleGaussianInteraction::computeDerivativeOfu(int particleNumber){

    int                     numberOfParticles       = m_system->getNumberOfParticles();
    int                     numberOfDimentions      = m_system->getNumberOfDimensions();
    double                  uDerivative             = 0;
    double                  a                       = m_system->getHardCoreDiameter();
    double                  uDoubleDerivative       = 0;
    double                  uTotalDoubleDerivative  = 0;
    std::vector <double>    uTotalDerivative        (numberOfDimentions);
    std::vector <double>    uAllStuff               (numberOfDimentions+1);

    auto m_particles = m_system->getParticles();
    auto ri         = m_particles[particleNumber]->getPosition();                              // (x_i, y_i, z_i)
    auto difference = getDistances(m_particles)[particleNumber];
    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        if (particleNumber != l1){
            auto rj         = m_particles[l1]->getPosition();                          // (x_j, y_j, z_j)
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

std::vector<double> SimpleGaussianInteraction::computeDerivativeOneParticle(int particleIndex){
    
    int                 numberOfDimensions          = m_system->getNumberOfDimensions();
    std::vector<double> derivativeVector            (numberOfDimensions);
    auto m_particles = m_system->getParticles();

    auto r = m_particles[particleIndex]->getPosition();

    for (int j3=0; j3<numberOfDimensions; j3++){
        if (j3 == 2){
            derivativeVector[j3] = -2*getParameters()[0]*m_beta*r[j3];
        }else{
            derivativeVector[j3] = -2*getParameters()[0]*r[j3];
        }
    }

    return derivativeVector;
    
}

bool SimpleGaussianInteraction::calculateInterparticleDistances(std::vector <class Particle*> particles){
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector <std::vector <double>>                distances       (numberOfParticles); // a matrix of the distance between all particles 
    std::vector <double>                              difference      (numberOfParticles); // the distance between particle j and all other particles i where j>i
    std::vector <double>                              vectorDistance  (numberOfDimensions); 

    double a = m_system->getHardCoreDiameter();
    
    std::vector <class Particle*> m_particles = particles;

    int distanceCheck = 0;
    std::vector <double> r1(numberOfDimensions), r2(numberOfDimensions);

    for (int k1 = 0; k1 < numberOfParticles; k1++){
        // std::cout << "entering loop. Particle 1 is " << k1 << std::endl;

        r1 = m_particles[k1]->getPosition();
        // std::cout << "got position" << std::endl;
        for (int k2 = 0; k2 <numberOfParticles; k2++){
            // std::cout << "entering loop. Particle 2 is " << k2 << std::endl;
            // std::cout << "hence calculating: r_" << k1 << k2 << std::endl;
            if (k1 != k2){
                r2 = m_particles[k2]->getPosition();
                for (int k3 = 0; k3<numberOfDimensions; k3++){
                    difference[k2] +=  (r1[k3]-r2[k3])*(r1[k3]-r2[k3]);              // (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2
                }
                difference[k2] = sqrt(difference[k2]);                                  // sqrt((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2)

                // std::cout << "this is r_ij for p " << k1+1 << std::endl;
                // std::cout << difference[0] << ",\t" << difference[1] << ",\t" << std::endl;
        

                if (difference[k2] < a){
                    // std::cout << "too small difference: " << difference[k2] << "<" << a << std::endl;
                    distanceCheck +=1;
                }
            }else{
                difference[k2] = 0;
            }
        }
        // std::cout << "all distances are big enough" << std::endl;

        distances[k1] = difference;

        // std::cout << "have saved the vector of distance to the matrix of distances ones" << std::endl;

    }
    setDistances(distances);

    if (distanceCheck > 0){
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

double SimpleGaussianInteraction::calculatePositionSumSquared(){

    // std::cout << "skal brukes" << std::endl;

    double vectorSumSquared = 0.0;
    auto m_particles = m_system->getParticles();
    
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    for (int i10 = 0; i10<numberOfParticles;i10++){
        auto r = m_particles[i10]->getPosition();
        for (int n10=0; n10<numberOfDimensions; n10++){
            if (n10 == 2){
                vectorSumSquared += m_beta*r[n10]*r[n10];
            }else{
                vectorSumSquared += r[n10]*r[n10];
            }
        }
    }

    // std::cout << "vector sum:" << vectorSumSquared << std::endl;

    return vectorSumSquared;
}