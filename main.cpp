#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

using namespace std;


int main() {

    double stopCriteria     = 1e-9;         // Stopping criteria for energy vs exact energy.
    double minimizationRate = 0.05;          // Model for double derivative of energy.
    double energyDerivative = 1.0;
    double alpha            = 0.45;          // Variational parameter - first guess.
    double alphaChange      = 0.5;
    double alphaNew = alpha;

    for (int k=0;  alphaChange > stopCriteria; k++){
        // clock_t start, end;
        // /* Recording the starting clock tick.*/
        // start = clock(); 

        std::cout << " alpha: " << alpha << endl;
        std::cout << " derivative: " << energyDerivative << endl;

        int numberOfDimensions  = 1;
        int numberOfParticles   = 1;
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double stepLength       = 1;          // Metropolis step length.
        double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
        double timeStep         = 1e-2;

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (true);
        system->setImportanceSampling       (true, timeStep);
        system->runMetropolisSteps          (numberOfSteps);
        
        energyDerivative = system->getHamiltonian()->computeEnergyDerivative(system->getParticles());

        alphaNew = alpha - minimizationRate*energyDerivative/numberOfParticles;

        alphaChange = std::abs(alphaNew-alpha);
        alpha = alphaNew;

    //     end = clock(); 
    
    //     // Calculating total time taken by the program. 
    //     double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    //     std::cout << "Time taken by program is : " << fixed  
    //         << time_taken; 
    //     std::cout << " sec " << std::endl; 
    }

    return 0;
}
