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

    // if (m_system->getSteps() == 4532 || m_system->getSteps() == 4531){
    //     std::cout << "positionSum: " << rSum << std::endl;
    // }

    double interactionPart = evaluateCorrelationPart();

    // if (m_system->getSteps() == 4532 || m_system->getSteps() == 4531){
    //     std::cout << "interactionPart: " << interactionPart << std::endl;
    // }

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
    

    for (int j1 = 0; j1 < numberOfParticles-1; j1++){

        for (int j2 = j1+1; j2 <numberOfParticles; j2++){
            distances_j1_j2 = distances[j1][j2];
            // std::cout << "distance: " << distances_j1_j2 << std::endl;
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

    // if ( m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "step: " << m_system->getSteps()-1e5 << std::endl;
    //     std::cout << "rSum2: " << rSum2 << std::endl;
    // }

    interactionPart = computeInteractionPartOfDoubleDerivative();

    // if (m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "step: " << m_system->getSteps()-1e5 << std::endl;
    //     std::cout << "Interaction part: " << interactionPart << std::endl;
    // }
    
    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2) + interactionPart;
}

std::vector<double> SimpleGaussianInteraction::computeDerivative(int particleIndex){
    
    int                     numberOfDimensions          = m_system->getNumberOfDimensions();
    double                  derivative_psi_ob           = 0;
    std::vector <double>    vectorWithInteraction       (numberOfDimensions);

    auto m_particles         = m_system->getParticles();
    auto ri                  = m_particles[particleIndex]->getPosition();
    auto derivative_psi_in   = computeDerivativeOfu(particleIndex);

    for (int n8=0; n8<numberOfDimensions; n8++){
        if (n8 == 2){
            derivative_psi_ob = -2*getParameters()[0]*ri[n8]*m_beta;
        }else{
            derivative_psi_ob = -2*getParameters()[0]*ri[n8];
        }
        vectorWithInteraction[n8] = derivative_psi_ob + derivative_psi_in[n8];
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
    int         numberOfDimensions = m_system->getNumberOfDimensions();
    auto        m_particles = m_system->getParticles();


    for(int l5 = 0; l5 < numberOfParticles; l5++){ 

        auto derivativePhi = computeDerivativeOneParticle(l5);

        auto uPart = computeDerivativeOfu(l5);

        for (int l4 = 0; l4<numberOfDimensions; l4++){
            firstTerm += derivativePhi[l4]*uPart[l4];
            secondTerm += uPart[l4]*uPart[l4];
        }
        thirdTerm += uPart[numberOfDimensions+1];
        
    }

    // if (m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "step: " << m_system->getSteps()-1e5 << std::endl;
    //     std::cout << "firstTerm: " << firstTerm << std::endl;
    //     std::cout << "secondTerm: " << secondTerm << std::endl;
    //     std::cout << "thirdTerm: " << thirdTerm << std::endl;
    // }

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
    auto difference = getDistances(m_particles);
    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        if (particleNumber != l1){
            auto rj         = m_particles[l1]->getPosition();                          // (x_j, y_j, z_j)
            auto rLength    = difference[particleNumber][l1];                                   // r_ij
    
        //     if (m_system->getSteps() == 4533+(int)1e5){
        //         std::cout << "particle: " << l1 << std::endl;
        //         std::cout << "rj: " << rj[0] <<", "<< rj[1]<<", " <<rj[2] << std::endl;
        //         std::cout << "rLength: " << rLength << std::endl;
        // }

            /* Here sum u'(r_ij) is determined based on the relationship 
            between r_ij (distance between particles) and a (hard core diameter) */

            if (rLength <= a){
                uDerivative = -1e50;
                uDoubleDerivative = -1e50;
                std::cout << "r_ij < a in computeDerivativeOfu" << std::endl;
            }else{
                uDerivative = -a/(a*rLength-rLength*rLength);
                uDoubleDerivative = a*(a-2*rLength)/(rLength*rLength*(a-rLength)*(a-rLength));

                // if (m_system->getSteps() == 4533+(int)1e5){
                //     std::cout << "uDerivative: " << uDerivative << std::endl;
                //     std::cout << "uDoubleDerivative: " << uDoubleDerivative << std::endl;
                // }
            }

            for (int l3 = 0; l3<numberOfDimentions; l3++){
                uTotalDerivative[l3] += ((ri[l3]-rj[l3])/rLength)*uDerivative;
            }

            uTotalDoubleDerivative += uDoubleDerivative + ((numberOfDimentions-1)/rLength)*uDerivative;
        }
    }

    for (int l4 = 0; l4<numberOfDimentions; l4++){
        uAllStuff[l4] = uTotalDerivative[l4];
        // if (m_system->getSteps() == 4533+(int)1e5){
        //     std::cout << "dim: " << l4 << std::endl;
        //     std::cout << "uTotalDerivative: " << uTotalDerivative[l4] << std::endl;
        // }
    }
    uAllStuff[numberOfDimentions+1] = uTotalDoubleDerivative;
    // if (m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "uTotalDoubleDerivative: " << uTotalDoubleDerivative << std::endl;
    // }

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

    std::vector <std::vector <double> > distances(numberOfParticles, std::vector<double>(numberOfParticles, (double) 0)); // a matrix of the distance between all particles 
    
    double a = m_system->getHardCoreDiameter();
    
    std::vector <class Particle*> m_particles = particles;

    int distanceCheck = 0;
    std::vector <double> r1(numberOfDimensions), r2(numberOfDimensions);

    int times = 0;

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
                    distances[k1][k2] +=  (r1[k3]-r2[k3])*(r1[k3]-r2[k3]);              // (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2
                }
                distances[k1][k2] = sqrt(distances[k1][k2]);                                  // sqrt((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2)
                // std::cout << "this is r_ij for p " << k1+1 << std::endl;
                // std::cout << difference[0] << ",\t" << difference[1] << ",\t" << std::endl;

                if (distances[k1][k2] < a){
                    // std::cout << "too small difference: " << difference[k2] << "<" << a << std::endl;
                    distanceCheck +=1;
                }
            }else{
                distances[k1][k2] = 0;
            }
        }
        // std::cout << "all distances are big enough" << std::endl;

        // std::cout << "have saved the vector of distance to the matrix of distances ones" << std::endl;

    }

    // if (m_system->getSteps() == 4533+(int)1e5){
    //     std::cout << "[ " << distances[0][0]<< " , " << distances[0][1] << " , " << distances[0][2] << " ] " << std::endl;
    //     std::cout << "[ " << distances[1][0]<< " , " << distances[1][1] << " , " << distances[1][2] << " ] " << std::endl;
    //     std::cout << "[ " << distances[2][0]<< " , " << distances[2][1] << " , " << distances[2][2] << " ] " << std::endl;
    // }

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