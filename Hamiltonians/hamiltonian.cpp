#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

/*
double Hamiltonian::computeDoubleDerivativeNumerically(std::vector<class Particle*> particles) {
    std::vector<Particle*> particles_more = particles;
    std::vector<Particle*> particles_less = particles;
    std::vector<Particle*> particles_org = particles;
    
    double step = 0.001;
    double dpsidr2 = 0.0;
    double rSum2 = 0.0;

    int numberOfParticles = particles.size();
    int numberOfDimensions = particles[0]->getPosition().size();

    std::vector<double> r(numberOfDimensions);

    // Finding the sum for r - original
    for(int i3=0; i3<numberOfParticles; i3++){
        r = particles[i3]->getPosition();
        for(int n3=0; n3<numberOfDimensions;n3++){
            rSum2 += r[n3]*r[n3];
        }
    }

    // Changing the copies to obtain an array of "next step" and "previous step"
    for(int i4=0; i4<numberOfParticles; i4++){
         for(int n4=0; n4<numberOfDimensions;n4++){
            // finding the next step and previous step for one particle in one spesific dimension
            particles_less[i4]->adjustPosition(step, n4);
            particles_more[i4]->adjustPosition(step, n4);

            dpsidr2 += (m_system->getWaveFunction()->evaluate(particles_more)+m_system->getWaveFunction()->evaluate(particles_less))/(step*step);

            // Resetting the arrays so that a new particle and spesific dimension can be calculated
            particles_less[i4]->adjustPosition(-step, n4);
            particles_more[i4]->adjustPosition(-step, n4);
        }
    }
    


    dpsidr2 += -2*numberOfParticles*numberOfDimensions*m_system->getWaveFunction()->evaluate(particles_more)/(step*step);

return dpsidr2


    
    
}*/
