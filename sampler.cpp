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
#include <omp.h>

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

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    
    // Sampling if the equilibrium stage is passed
    double localEnergy = 0;

    if (m_stepNumber > m_system->getEquilibration()){
        if (acceptedStep == 1) { m_numberOfAcceptedSteps += 1; }
        localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
        m_cumulativeEnergy  += localEnergy;
        m_cumulativeEnergySquared += localEnergy*localEnergy;
        m_cumulativeEnergyDerivative += localEnergy*m_system->getWaveFunction()->computeAlphaDerivative();
        m_cumulativeAlphaDerivative += m_system->getWaveFunction()->computeAlphaDerivative();
    }

    m_stepNumber++;
}

void Sampler::sampleAllEnergies(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeEnergyDerivative = 0;
        m_cumulativeAlphaDerivative = 0;
        m_localEnergyVector = std::vector <double>();
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    
    // Sampling if the equilibrium stage is passed

        if (m_stepNumber >= m_system->getEquilibration()){
            if (acceptedStep == 1) { m_numberOfAcceptedSteps += 1; }

            m_system->getWaveFunction()->updateOneBodyDensity();
            double localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
            
            m_cumulativeEnergy  += localEnergy;
            m_localEnergyVector.push_back(localEnergy);
            m_cumulativeEnergySquared += localEnergy*localEnergy;
            m_cumulativeEnergyDerivative += localEnergy*m_system->getWaveFunction()->
                                            computeAlphaDerivative();
            m_cumulativeAlphaDerivative += m_system->getWaveFunction()->computeAlphaDerivative();
           
        }

        m_stepNumber++;
    
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibration();
    double  ac = 100*m_numberOfAcceptedSteps/ms;
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(ef) << endl;
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

void Sampler::printOutputToEnergyAlphaFile(){
    
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    double  ef = m_system->getEquilibration();

    ofstream myfile;
    double alpha = m_system->getWaveFunction()->getParameters()[0];

    string filename = "Output/" + m_system->getFileName() + "_energy_alpha.txt";

    if (m_firstCriteria == 0) { 
        
        myfile.open (filename, ios::out | ios::trunc);
        myfile << "  -- System info -- " << endl;
        myfile << " Number of particles  : " << np << endl;
        myfile << " Number of dimensions : " << nd << endl;
        myfile << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
        myfile << " Number of equilibration steps  : 10^" << std::log10(ef) << endl;
        myfile << " Step length/ time step: " << m_system->getStepLength() << endl;
        myfile << "===================================================" << endl;
        myfile << "Energy: \t Alpha: \t Acceptance: \n"; 
        myfile << m_energy << "\t" << alpha << "\t" << 100.0*(double)m_numberOfAcceptedSteps/(double)m_numberOfMetropolisSteps <<"\n";
        myfile.close(); 
    }else{
        // std::cout << "printing" << std::endl;
    myfile.open (filename, ios::out | ios::app);
        
    myfile << m_energy << "\t" << alpha << "\t" << 100.0*(double)m_numberOfAcceptedSteps/(double)m_numberOfMetropolisSteps <<"\n";
    myfile.close();
    std::cout << m_energy << "\t" << alpha <<"\n";
    }
}

void Sampler::printOutputToEnergyFile(){

    ofstream myfile;

    evaluateNumberOfCyclesIncluded();
    std::cout << "accepted steps: " << (double)m_numberOfAcceptedSteps/(double)m_numberOfCyclesIncluded << endl;

    string filename = m_system->getFileName() + "_energy.txt";

    // cout << "energy vector length: " << m_localEnergyVector.size() << endl;
    // cout << "should be: " << (int) pow(2,20) << endl;

    myfile.open (filename, ios::out | ios::trunc);

    myfile << " # mean: " << m_energy << "\n";
    myfile << " # std: " << "sqrt(" << m_energySquared << " - " << m_energy << "^2) = " << sqrt(m_energySquared - m_energy*m_energy) << "\n";


    for (int n3 = 0; n3<(int) m_localEnergyVector.size(); n3++){
        myfile << m_localEnergyVector[n3] << "\n";
    }
    myfile.close();
}

void Sampler::printOneBodyDensityToFile(){
    
    std::vector <std::vector <double>> oneBodyDensity = m_system->getWaveFunction()->getOneBodyDensity();
    evaluateNumberOfCyclesIncluded();

    ofstream myfile;

    string filename = m_system->getFileName() + "_density.txt";

    myfile.open (filename, ios::out | ios::trunc);

    int numberOfParticles = (double) m_system->getNumberOfParticles();

    for (int n5 = 0; n5 < (int)oneBodyDensity[0].size(); n5++){
        myfile << oneBodyDensity[0][n5] << "\t";
        for (int m5 = 1; m5 < m_system->getNumberOfDimensions()+1; m5++){
            myfile << oneBodyDensity[m5][n5]/((double) m_numberOfCyclesIncluded*(double) numberOfParticles) << "\t";
        }
        myfile << "\n";
    }

    myfile.close();
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     * Take away the non-physical stuff before eqilibrium - so not all steps
     */

    evaluateNumberOfCyclesIncluded();

    cout << "number of cycles: " << m_numberOfCyclesIncluded << endl;

    m_energy = (m_cumulativeEnergy / (double) m_numberOfCyclesIncluded);
    m_energySquared = m_cumulativeEnergySquared / (double) m_numberOfCyclesIncluded;

    cout << "energySquared = " << m_cumulativeEnergySquared << " / " << m_numberOfCyclesIncluded << endl;

    m_derivative = 2* m_cumulativeEnergyDerivative / (double) m_numberOfCyclesIncluded
                    - 2*(m_cumulativeAlphaDerivative / (double) m_numberOfCyclesIncluded)*m_energy;
}

void Sampler::setFileOutput(int firstCriteria){
    m_firstCriteria = firstCriteria;
    // m_filename = filename;
}

void Sampler::evaluateNumberOfCyclesIncluded(){
    m_numberOfCyclesIncluded = m_numberOfMetropolisSteps - m_system->getEquilibration();
}
