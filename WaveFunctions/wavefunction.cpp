#include "wavefunction.h"
#include "system.h"
#include "particle.h"
#include <cmath>


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setOneBodyDensityBins(int numberOfBins, double densityLength){
    m_densityLength = densityLength;
    m_numberOfBins = numberOfBins;

    for (int n3 = 0; n3<numberOfBins; n3++){
        m_oneBodyDensity.push_back(0);
    }
}

void WaveFunction::updateOneBodyDensity(std::vector<class Particle*> particles){
    
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double binLength = m_densityLength/(double) m_numberOfBins;

    for (int l = 0; l < numberOfParticles; l++){
        double rLength = 0;
        auto r = particles[l]->getPosition();
        for (int j3 = 0; j3 < numberOfDimensions; j3++){
            if (j3 == 2){
                rLength += r[j3]*r[j3];
            }
        }
        rLength = sqrt(rLength);
        for (int j5 = 0; j5 < m_numberOfBins; j5++){
            if (rLength >= j5*binLength && rLength < (j5+1)*binLength){
                m_oneBodyDensity[j5] += 1;
            }
        }
    }


}
