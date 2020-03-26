#include "wavefunction.h"
#include "system.h"
#include "particle.h"
#include <cmath>
#include <iostream>


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}

void WaveFunction::setOneBodyDensityBins(int numberOfBins, double densityLength){
    m_densityLength = densityLength;

    if (numberOfBins % 2 == 0){
        m_numberOfBins = numberOfBins;
    }else{
        std::cout << "Number of bins must be able to be divided by 2" << std::endl;
    }
    
    int numberOfDimensions = m_system->getNumberOfDimensions();
    std::vector<double> bins(m_numberOfBins);
    std::vector<double> binStartValue(m_numberOfBins);

    for (int n3 = 0; n3<numberOfBins; n3++){
        bins[n3] = 0;
        binStartValue[n3] = (-m_numberOfBins/2+n3)*(densityLength/(double)numberOfBins);
    }
    // Creating an array to have control over what values the bins represent
    m_oneBodyDensity.push_back(binStartValue);

    for (int n6 = 0; n6<numberOfDimensions; n6++){
        m_oneBodyDensity.push_back(bins);
    }

}

void WaveFunction::updateOneBodyDensity(std::vector<class Particle*> particles){
    
    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double binLength = m_densityLength/(double) m_numberOfBins;

    // std::cout << "number of bins: " << m_numberOfBins << std::endl;
    // std::cout << "length of bins: " << m_oneBodyDensity[0].size() << std::endl;

    // std::cout << "check bin place: " << m_oneBodyDensity[0][302] << std::endl;


    for (int l = 0; l < numberOfParticles; l++){

        auto r = particles[l]->getPosition();

        for (int j3 = 0; j3 < numberOfDimensions; j3++){
                // std::cout << "r is: " << r[j3] << std::endl;
            for (int j5 = -m_numberOfBins/2; j5 < m_numberOfBins/2; j5++){
                // std::cout << j5 << std::endl;
            
                if (j5 < 0){
                    if (r[j3] > (j5)*binLength && r[j3] <= (j5+1)*binLength){
                        // std::cout << "r_min: " << (j5)*binLength << " r_max: " << (j5+1)*binLength << std::endl;
                        m_oneBodyDensity[j3+1][j5+m_numberOfBins/2] += 1;
                    }
                }else{
                    if (r[j3] >= (j5)*binLength && r[j3] < (j5+1)*binLength){
                        // std::cout << "r_min: " << (j5)*binLength << " r_max: " << (j5+1)*binLength << std::endl;
                        m_oneBodyDensity[j3+1][j5+m_numberOfBins/2] += 1;
                    }
                }
                // std::cout << "this is ok. d = "<< j3 << " bin =" << j5+m_numberOfBins/2 << std::endl;

            }
        }
    }
}
