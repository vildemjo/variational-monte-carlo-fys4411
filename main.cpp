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

    std::vector<int> Ns = {500};

    for (int N = 0; N < Ns[-1]; N++){
        clock_t start, end;
        /* Recording the starting clock tick.*/
        start = clock(); 

        int numberOfDimensions  = 1;
        int numberOfParticles   = Ns[N];
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double alpha            = 0.5;          // Variational parameter.
        double stepLength       = 0.1;          // Metropolis step length.
        double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
        double timeStep         = 1e-2;

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (false);
        system->setImportanceSampling       (false, timeStep);
        system->runMetropolisSteps          (numberOfSteps);
        
        end = clock(); 
    
        // Calculating total time taken by the program. 
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
        std::cout << "Time taken by program is : " << fixed  
            << time_taken; 
        std::cout << " sec " << std::endl; 
        
    }
    return 0;
}
