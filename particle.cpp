#include "particle.h"
#include <cassert>
#include "system.h"
#include "InitialStates/initialstate.h"

Particle::Particle(System* system) {
    m_system = system;
}

void Particle::setPosition(const std::vector<double> &position) {
    assert(position.size() == m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
    if (m_system->getInteractionOrNot() == true){
        m_system->getInitialState()->calculateInterparticleDistances();
    }
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
