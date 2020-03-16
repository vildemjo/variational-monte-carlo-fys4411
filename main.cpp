#include <iostream>
#include <fstream>
#include "system.h"
#include "particle.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <string>

using namespace std;

void gradientDecentRun(string filename, double alphaStart, double minimizationRate, double stopCriteria, bool analytic, bool importance, bool interactionOrNot, double hardCoreDiameter);
void alphaListRun(string filename, double alphaStart, double alphaStop, double alphaStep, bool analytic, bool importance, bool interactionOrNot, double hardCoreDiameter);

int main() {

    clock_t start, end;
    /* Recording the starting clock tick.*/
    start = clock();

    bool analyticOrNot = true;
    bool importanceOrNot = false;
    bool interactionOrNot = false;
    double hardCoreDiameter = 0.0043;

    bool gradientDescent = false;
    double minimizationRate = 0.1;
    double alphaGuess = 0.4;
    double stopCriteria = 1e-9;

    bool alphaList = true;
    double alphaStart = 0.3;
    double alphaStop = 0.01;
    double alphaStep = 0.02;
    string filename = "with_int/num_2p_1d_a_0-0043_energy_alpha.txt";

    // elliptical or not
    double beta = 2.82843;    

    if (gradientDescent == true){
        gradientDecentRun(filename, alphaGuess, minimizationRate, stopCriteria, analyticOrNot, importanceOrNot, interactionOrNot, hardCoreDiameter);
    }
    if(alphaList == true){
        alphaListRun(filename, alphaStart, alphaStop, alphaStep, analyticOrNot, importanceOrNot, interactionOrNot, hardCoreDiameter);
    }

    end = clock(); 

    // Calculating total time taken by the program. 
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    std::cout << "Time taken by program is : " << fixed  
        << time_taken; 
    std::cout << " sec " << std::endl; 

    return 0;
}

void gradientDecentRun(string filename, 
                    double alphaGuess, 
                    double minimizationRate, 
                    double stopCriteria, 
                    bool analyticOrNot, 
                    bool importanceOrNot, 
                    bool interactionOrNot, 
                    double hardCoreDiameter){

    int firstCriteria               = 0;                  // Criteria to print header in outputfile

    double alpha                    = alphaGuess;         // Variational parameter first guess
    double m_stopCriteria           = stopCriteria;       // Stopping criteria for energy vs exact energy.
    double m_minimizationRate       = minimizationRate;   // Model for double derivative of energy.
    double energyDerivative         = 1.0;
    double alphaChange              = 0.5;
    double energyChange             = 1;
    double alphaNew                 = alpha;
    double energyNew                = 0.8; 
    double gradientPrintToFileOrNot = false;

    ofstream file;
    if (gradientPrintToFileOrNot == true){
        file.open ("Output/gradient_descent.txt", ios::out | ios::trunc);
        file << "Alpha: \t Energy: \t Derivative: \n";
        file.close();
    }

    for (int k=0;  energyChange > m_stopCriteria; k++){

        int numberOfDimensions  = 1;
        int numberOfParticles   = 1;
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double stepLength       = 0.5;          // Metropolis step length.
        double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
        double timeStep         = 1e-2;
        bool printToFileOrNot   = false;
        double energy           = 0;
    

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInteractionOrNot         (interactionOrNot, hardCoreDiameter); // Must come before initialize state
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (analyticOrNot);
        system->setImportanceSampling       (importanceOrNot, timeStep);
        system->runMetropolisSteps          (numberOfSteps, printToFileOrNot, firstCriteria);
    
        firstCriteria = 1;
        
        energyNew = system->getSampler()->getEnergy();
        energyDerivative = system->getSampler()->getDerivative();
        alphaNew = alpha - m_minimizationRate*energyDerivative/numberOfParticles;
        
        energyChange = std::abs(energyNew - energy);
        alpha = alphaNew;
        energy = energyNew;

        if (gradientPrintToFileOrNot == true){
            file.open ("Output/gradient_descent.txt", ios::out | ios::app);
            file << alpha << "\t" << energy << "\t" << energyDerivative << "\n";
            file.close();
        }
        
    }
    

}

void alphaListRun(string filename,
                double alphaStart,
                double alphaStop, 
                double alphaStep, 
                bool analytic, 
                bool importanceOrNot,
                bool interactionOrNot, 
                double hardCoreDiameter){

    int firstCriteria = 0;                  // Criteria to print header in outputfile

    double alpha = alphaStart;
    for(int a=0; alpha > alphaStop+alphaStep; a++){

        double stopCriteria     = 1e-9;         // Stopping criteria for energy vs exact energy.
        int numberOfDimensions  = 1;
        int numberOfParticles   = 10;
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double stepLength       = 0.5;          // Metropolis step length.
        double equilibration    = 0.01;          // Fraction of the total steps used for equilibration
        double timeStep         = 1e-2;
        bool printToFileOrNot   = true;

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInteractionOrNot         (interactionOrNot, hardCoreDiameter);
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (analytic);
        system->setImportanceSampling       (importanceOrNot, timeStep);
        system->setFileName                 (filename);
        system->runMetropolisSteps          (numberOfSteps, printToFileOrNot, firstCriteria);
            
        alpha -= alphaStep;
        firstCriteria = 1;
    }
 
}