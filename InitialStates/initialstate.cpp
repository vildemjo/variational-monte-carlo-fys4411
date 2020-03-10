#include "initialstate.h"
#include <cmath>
#include "../particle.h"
#include "../system.h"
#include <cassert>
#include "../Hamiltonians/hamiltonian.h"

InitialState::InitialState(System* system) {
    m_system = system;
}

void InitialState::calculateInterparticleDistances(){
    std::vector <std::vector <double>>      distances  (m_numberOfParticles);
    std::vector <double>                    difference (m_numberOfParticles);

    double a = m_system->getHardCoreDiameter();
    m_system->getHamiltonian()->setInteractionPotential(false);

    for (int j1 = 0; j1 < m_numberOfParticles-1; j1++){
        auto r1 = m_particles[j1]->getPosition();
        for (int j2 = j1+1; j2 <m_numberOfParticles; j2++){
            auto r2 = m_particles[j2]->getPosition();
            for (int j3 = 0; j3<m_numberOfDimensions; j3++){
                difference[j2] +=  r1[j3]-r2[j3];
            }
            difference[j2] = sqrt(difference[j2]);
            if (difference[j2] < a){
                m_system->getHamiltonian()->setInteractionPotential(true);  // Telling the interaction potential that a distance is smaller than a
            }
        }
        distances[j1] = difference;
    }

    setDistances(distances);
}

void InitialState::setDistances(std::vector<std::vector<double>> distances){
    m_distances = distances;
}

void InitialState::updateDistances(int particleNumber){
    assert(particleNumber < m_numberOfParticles);
    
    std::vector <double> difference(m_numberOfParticles);
    auto r1 = m_particles[particleNumber]->getPosition();
    double a = m_system->getHardCoreDiameter();

    for(int j4 = particleNumber+1; j4<m_numberOfParticles; j4++){
        auto r2 = m_particles[particleNumber]->getPosition();
        for (int j5 = 0; j5<m_numberOfDimensions; j5++){
            difference[j4] += r1[j5]-r2[j5];
        }
        difference[j4] = sqrt(difference[j4]);
        if (difference[j4] < a){
            m_system->getHamiltonian()->setInteractionPotential(true);  // Telling the interaction potential that a distance is smaller than a
        }
    }
    m_distances[particleNumber] = difference;
}