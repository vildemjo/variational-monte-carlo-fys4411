#include "initialstate.h"
#include <cmath>
#include "../particle.h"
#include "../system.h"
#include <cassert>
#include "../Hamiltonians/hamiltonian.h"
#include <iostream>

InitialState::InitialState(System* system) {
    m_system = system;
}

bool InitialState::calculateInterparticleDistances(){
    std::vector <std::vector <double>>                distances       (m_numberOfParticles); // a matrix of the distance between all particles 
    std::vector <double>                              difference      (m_numberOfParticles); // the distance between particle j and all other particles i where j>i
    std::vector <double>                              vectorDistance  (m_numberOfDimensions); 

    double a = m_system->getHardCoreDiameter();
    m_system->getHamiltonian()->setInteractionPotential(false);

    for (int j1 = 0; j1 < m_numberOfParticles; j1++){
        // std::cout << "entering loop. Particle 1 is " << j1 << std::endl;

        auto r1 = m_particles[j1]->getPosition();
        for (int j2 = 0; j2 <m_numberOfParticles; j2++){
            // std::cout << "entering loop. Particle 2 is " << j2 << std::endl;
            // std::cout << "hence calculating: r_" << j1 << j2 << std::endl;
            if (j1 != j2){
                auto r2 = m_particles[j2]->getPosition();
                for (int j3 = 0; j3<m_numberOfDimensions; j3++){
                    difference[j2] +=  (r1[j3]-r2[j3])*(r1[j3]-r2[j3]);              // (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2
                }
                difference[j2] = sqrt(difference[j2]);                                  // sqrt((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2)

                // std::cout << "this is r_ij for p " << j1+1 << std::endl;
                // std::cout << difference[0] << ",\t" << difference[1] << ",\t" << difference[2] << ",\t" << std::endl;
        

                if (difference[j2] < a){
                    // std::cout << "too small difference" << std::endl;
                    m_system->getHamiltonian()->setInteractionPotential(true);  // Telling the interaction potential that a distance is smaller than a
                    return false;
                }
            }
        }
        // std::cout << "all distances are big enough" << std::endl;

        distances[j1] = difference;

        // std::cout << "have saved the vector of distance to the matrix of distances ones" << std::endl;

    }
    setDistances(distances);
    // std::cout << "have set the distances" << std::endl;
    return true;
}


void InitialState::setDistances(std::vector<std::vector<double>> distances){
    m_distances = distances;
    // std::cout << "works inside setDistance function too" << std::endl;
}


void InitialState::updateDistances(int particleNumber){
    
    bool checkDistance = calculateInterparticleDistances();
    if (checkDistance == false){
        // std::cout << "the distance got to small during stepping" << std::endl;
    }
    // assert(particleNumber < m_numberOfParticles);
    
    // std::vector <double> difference(m_numberOfParticles);

    // auto r1 = m_particles[particleNumber]->getPosition();
    // double a = m_system->getHardCoreDiameter();
    // int checkDistance = 0;

    // for(int j4 = 0; j4<m_numberOfParticles; j4++){
    //     if (j4 != particleNumber){
    //         auto r2 = m_particles[particleNumber]->getPosition();
    //         for (int j5 = 0; j5<m_numberOfDimensions; j5++){
    //             difference[j4] += r1[j5]-r2[j5];
    //         }
    //         difference[j4] = sqrt(difference[j4]);
    //         if (difference[j4] < a){
    //             checkDistance += 1;  // Checking wether one of the distances are smaller than a
    //         }
    //     }
    // }
    // m_distances[particleNumber] = difference;
    // if (checkDistance > 0){
    //     m_system->getHamiltonian()->setInteractionPotential(true);
    // } // If not the interaction potential stays as it was after initialization
}

std::vector<double> InitialState::evaluateDifferenceVector(){

    std::vector <double> differenceVector(m_numberOfDimensions);

    for (int k1 = 0; k1 < m_numberOfParticles-1; k1++){
        auto r1 = m_particles[k1]->getPosition();
        for (int k2 = k1+1; k2 <m_numberOfParticles; k2++){
            auto r2 = m_particles[k2]->getPosition();
            for (int k3 = 0; k3<m_numberOfDimensions; k3++){
                differenceVector[k3] += r1[k3]-r2[k3];
            }
        }
    }
    return differenceVector;
}