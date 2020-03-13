#include "initialstate.h"
#include <cmath>
#include "../particle.h"
#include "../system.h"
#include <cassert>
#include "../Hamiltonians/hamiltonian.h"

InitialState::InitialState(System* system) {
    m_system = system;
}

bool InitialState::calculateInterparticleDistances(){
    std::vector <std::vector <double>>                distances       (m_numberOfParticles); // a matrix of the distance between all particles 
    std::vector <double>                              difference      (m_numberOfParticles); // the distance between particle j and all other particles i where j>i
    std::vector <double>                              vectorDistance  (m_numberOfDimensions); 

    double a = m_system->getHardCoreDiameter();
    m_system->getHamiltonian()->setInteractionPotential(false);

    for (int j1 = 0; j1 < m_numberOfParticles-1; j1++){
        auto r1 = m_particles[j1]->getPosition();
        for (int j2 = j1+1; j2 <m_numberOfParticles; j2++){
            auto r2 = m_particles[j2]->getPosition();
            for (int j3 = 0; j3<m_numberOfDimensions; j3++){
                difference[j2] +=  (r1[j3]-r2[j3])*(r1[j3]-r2[j3]);              // (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2
            }
            difference[j2] = sqrt(difference[j2]);                                  // sqrt((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2)

            if (difference[j2] < a){
                m_system->getHamiltonian()->setInteractionPotential(true);  // Telling the interaction potential that a distance is smaller than a
                return false;
            }
        }
        distances[j1] = difference;
    }
    setDistances(distances);
    return true;
}


void InitialState::setDistances(std::vector<std::vector<double>> distances){
    m_distances = distances;
}

void InitialState::updateDistances(int particleNumber){
    assert(particleNumber < m_numberOfParticles);
    
    std::vector <double> difference(m_numberOfParticles);

    auto r1 = m_particles[particleNumber]->getPosition();
    double a = m_system->getHardCoreDiameter();
    int checkDistance = 0;

    for(int j4 = particleNumber+1; j4<m_numberOfParticles; j4++){
        auto r2 = m_particles[particleNumber]->getPosition();
        for (int j5 = 0; j5<m_numberOfDimensions; j5++){
            difference[j4] += r1[j5]-r2[j5];
        }
        difference[j4] = sqrt(difference[j4]);
        if (difference[j4] < a){
            checkDistance += 1;  // Telling the interaction potential that a distance is smaller than a
        }
    }
    m_distances[particleNumber] = difference;
    if (checkDistance > 0){
        m_system->getHamiltonian()->setInteractionPotential(true);
    } // If not it stays as it was after initialization (Should I add something that makes sure that it is ok after initialization?)
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