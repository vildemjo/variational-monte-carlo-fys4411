#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <string>

using namespace std;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
    }
    if (acceptedStep == 1) { m_numberOfAcceptedSteps += 1; }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    // Sampling if the equilibrium stage is passed
    double localEnergy;

    if (m_stepNumber > m_system->getEquilibrationFraction()*m_system->getNumberOfMetropolisSteps()){
        localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
        m_cumulativeEnergy  += localEnergy;
        m_cumulativeEnergySquared += localEnergy*localEnergy;
        m_cumulativeEnergyDerivative += localEnergy*m_system->getWaveFunction()->computeAlphaDerivative(m_system->getParticles());
        //m_cumulativeEnergyDerivative += m_system->getHamiltonian()->computeEnergyDerivative(m_system->getParticles());
        m_cumulativeAlphaDerivative += m_system->getWaveFunction()->computeAlphaDerivative(m_system->getParticles());

    }

    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    double  ac = 100*m_numberOfAcceptedSteps/ms;
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << " Accepted steps: " << ac << " %" << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Variance : " <<  m_energySquared - m_energy*m_energy << endl;
    cout << endl;
}

void Sampler::printOutputToFile(){
    
    ofstream myfile;
    double alpha = m_system->getWaveFunction()->getParameters()[0];

    string filename = "Output/" + m_system->getFileName();

    if (m_firstCriteria == 0) { 
        // Setting the correct filename for the different settings
        
        myfile.open (filename, ios::out | ios::trunc);
        myfile << "Energy: \t Alpha: \n"; 
        myfile.close(); 
    }
        
    myfile.open (filename, ios::out | ios::app);
        
    myfile << m_energy << "\t" << alpha << "\n";
    myfile.close();
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     * Take away the non-physical stuff before eqilibrium - so not all steps
     */

    int numberOfCyclesIncluded = (m_system->getNumberOfMetropolisSteps()*
                                    (1-m_system->getEquilibrationFraction()));

    m_energy = m_cumulativeEnergy / numberOfCyclesIncluded;
    m_energySquared = m_cumulativeEnergySquared / numberOfCyclesIncluded;
    m_derivative = 2*m_cumulativeEnergyDerivative / numberOfCyclesIncluded
                    - 2*(m_cumulativeAlphaDerivative / numberOfCyclesIncluded)*m_energy;
}

void Sampler::setFileOutput(int firstCriteria){
    m_firstCriteria = firstCriteria;
    // m_filename = filename;
}
