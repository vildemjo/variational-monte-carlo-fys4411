#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
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

    // std::cout << "I am inside evaluate wavefunction" << std::endl;

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    // std::cout << "there are " << numberOfParticles << " particles in the system" << std::endl;

    for(int i1=0; i1<numberOfParticles; i1++){
        auto r = particles[i1]->getPosition();
        // std::cout << "I can get the position of particle " << i1+1 << std::endl;
        for(int n1=0; n1<numberOfDimensions; n1++){
            rSum += r[n1]*r[n1];
        }
    }

    double interactionPart = 1;

    // std::cout << "test_ok" << std::endl;
   
    if (m_system->getInteractionOrNot() == true){

        auto a = m_system->getHardCoreDiameter();
        double uSum = 0;
        int uSumCheck = 0;
        std::vector<double> distances_j1(numberOfParticles);
        std::vector<double> distances_j2(numberOfParticles);
        double distances_j1_j2;

        auto distances = m_system->getInitialState()->getDistances();


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
        // Here the interaction part is set to zero directly (instead of exp(-infty))¨
        // if one of the distances are less than a
        if (uSumCheck > 0){
            interactionPart = 0;
        }else{
            interactionPart = exp(uSum);
        }
   }

    return exp(-m_parameters[0]*rSum)*interactionPart;

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

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double rSum2 = 0.0;
    double interactionPart = 0;

    for(int i2=0; i2<numberOfParticles; i2++){
        auto r = particles[i2]->getPosition();
        for(int n2=0; n2<numberOfDimensions; n2++){
            rSum2 += r[n2]*r[n2];
        }   
    }
    
    // std::cout << "The double derivative without interaction is okay" << std::endl;

    if (m_system->getInteractionOrNot() == true){
        interactionPart = computeInteractionPartOfDoubleDerivative(particles);
    }

    return (-2*m_parameters[0]*numberOfParticles*numberOfDimensions + 4*m_parameters[0]*m_parameters[0]*rSum2) + interactionPart;
}


std::vector<double> SimpleGaussian::computeDerivative(std::vector<class Particle*> particles){
    
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector<double> vectorSum(numberOfDimensions);

    if (m_system->getInteractionOrNot() == true){

        std::vector <double> vectorWithInteraction(numberOfDimensions);
        
        for (int i8 = 0; i8<m_system->getNumberOfParticles();i8++){

            // std::cout << "getting into loop works for p " << i8+1 << std::endl;

            auto uDerivative = computeDerivativeOfu(particles, i8);

            // std::cout << "computeDerivativeOfu works p is " << i8+1 << std::endl;

            auto r = particles[i8]->getPosition();

            // std::cout << "getting pos works for p " << i8+1 << std::endl;

            for (int n8=0; n8<numberOfDimensions; n8++){
                double phiPart = -2*getParameters()[0]*r[n8];
            
                // std::cout << "getting phiPart works for dim " << n8 << std::endl;

                vectorSum[n8] += phiPart;
                vectorWithInteraction[n8] += phiPart + uDerivative[n8];

                // std::cout << "getting uDerivative works for dim " << n8 << std::endl;
            }
        }

        return vectorWithInteraction;
    
    }else{
        for (int i8 = 0; i8<m_system->getNumberOfParticles();i8++){

            // std::cout << "getting into loop works for p " << i8+1 << std::endl;

            auto r = particles[i8]->getPosition();

            // std::cout << "getting pos works for p " << i8+1 << std::endl;

            for (int n8=0; n8<numberOfDimensions; n8++){
            double phiPart = -2*getParameters()[0]*r[n8];
            
            // std::cout << "getting phiPart works for dim " << n8 << std::endl;
            vectorSum[n8] += phiPart;
            }
        }
        return vectorSum;
    }

    

    
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


double SimpleGaussian::computeInteractionPartOfDoubleDerivative(std::vector<class Particle*> particles){

    double firstTerm; double secondTerm; double thirdTerm;
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimentions = m_system->getNumberOfDimensions();

    double uDerivative = 1;
    double uDoubleDerivative = 1;


    // Here we are summing over all the particles, to get the total 
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


std::vector <double> SimpleGaussian::computeDerivativeOfu(std::vector<class Particle*> particles, int particleNumber){
    
    // std::cout << "calulating derivative of u for p " << particleNumber+1 << std::endl;

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimentions = m_system->getNumberOfDimensions();

    double uDerivative = 0;
    double a = m_system->getHardCoreDiameter();

    std::vector <double> uTotalDerivative(numberOfDimentions);
    double uDoubleDerivative = 0;
    double uTotalDoubleDerivative = 0;

    std::vector <double> uAllStuff(numberOfDimentions+1);

    auto ri = particles[particleNumber]->getPosition();                              // (x_i, y_i, z_i)
    // std::cout << "getting position for p " << particleNumber+1 << std::endl;
    auto difference = m_system->getInitialState()->getDistances()[particleNumber];
    // std::cout << "getting distances for p " << particleNumber+1 << std::endl;

 
    // std::cout << "this is r_ij for p " << particleNumber+1 << std::endl;
    // std::cout << difference[0] << ",\t" << difference[1] << ",\t" << difference[2] << ",\t" << std::endl;
    

    for (int l1 = 0; l1 < numberOfParticles; l1++){
        // std::cout << "into loop for r_" << particleNumber << l1 << std::endl;
        if (particleNumber != l1){
            auto rj = particles[l1]->getPosition();                          // (x_j, y_j, z_j)
            // std::cout << "getting pos for r_" << particleNumber << l1 << std::endl;
            auto rLength = difference[l1];                                   // r_ij
            // std::cout << "getting |r_" << particleNumber << l1 << "|" << std::endl;

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

            uTotalDoubleDerivative += uDoubleDerivative + (2/rLength)*uDerivative;

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

std::vector<double> SimpleGaussian::computeDerivativeOneParticle(std::vector<class Particle*> particles, int particleIndex){
    
    int numberOfDimensions = m_system->getNumberOfDimensions();

    std::vector<double> derivativeVector(numberOfDimensions);

    auto r = particles[particleIndex]->getPosition();

    for (int n3=0; n3<numberOfDimensions; n3++){
        derivativeVector[n3] = -2*getParameters()[0]*r[n3];
    }

    return derivativeVector;
    
}