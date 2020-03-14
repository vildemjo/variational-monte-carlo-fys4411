#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../particle.h"
#include "../system.h"

using std::cout;
using std::endl;

RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
    // std::cout << "Done setting up in the class" << std::endl;
}

void RandomUniform::setupInitialState() {

    bool positionCheck = false;
    while ( positionCheck == false){

        for (int m3=0; m3 < m_numberOfParticles; m3++) {
            std::vector<double> position = std::vector<double>();
        
            for (int j=0; j < m_numberOfDimensions; j++) {

                position.push_back(m_system->getStepLength()*(Random::nextDouble()-0.5));
            }

            m_particles.push_back(new Particle(m_system));
            m_particles.at(m3)->setNumberOfDimensions(m_numberOfDimensions);
            m_particles.at(m3)->setPosition(position);
            m_particles.at(m3)->setParticleIndex(m3);
        }

        if (m_system->getInteractionOrNot() == true){
            positionCheck = calculateInterparticleDistances();
        }else{
            positionCheck = true;
        }
    }
    // std::cout << "Im out of the while-loop" << std::endl;
}
