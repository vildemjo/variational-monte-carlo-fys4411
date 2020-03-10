#include "particle.h"
#include <cassert>
#include "system.h"

Particle::Particle(System* system) {
    m_system = system;
}

void Particle::setPosition(const std::vector<double> &position) {
    assert(position.size() == m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position.at(dimension) += change;
    m_system->getInitialState()->updateDistances(m_particleIndex);
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
