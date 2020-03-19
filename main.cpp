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
#include <omp.h>

using namespace std;

void gradientDecentRun(string filename, double alphaStart, double minimizationRate, double stopCriteria, bool analytic, bool importance, bool interactionOrNot, double hardCoreDiameter, double beta);
void alphaListRun(string filename, double alphaStart, double alphaStop, double alphaStep, bool analytic, bool importance, bool interactionOrNot, double hardCoreDiameter, double beta);
string setMethodName(bool analyticOrNot);

int main() {

#pragma omp

/* The standard set-up */

    bool analyticOrNot      = true;
    double hardCoreDiameter = 0.0043;
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
    double timeStep         = 1e-2;
    bool printToFileOrNot   = true;
    int firstCriteria       = 0;            // print header in file
    int alpha               = 0.5;

    // elliptical or spherical trap (2.82843 or 1.0)
    double beta = 1.0;//2.82843;    // omega_normal^2/omega_ho^2

/* Exercise b */

// comparing with exact answear
/*

    double alphaStart = 0.7;
    double alphaStop = 0.01;
    double alphaStep = 0.02;

    int firstCriteria = 0;                  // Criteria to print header in outputfile
    bool analyticOrNot = true;
    double beta = 1;

    string methodName = setMethodName(analyticOrNot);

    std::vector <int> Ns = {1, 10, 100, 500};
    std::vector <int> ds = {1, 2, 3};

    for (int d = 0; d < ds.size(); d++){

        for (int N = 0; N < Ns.size(); N++ ){

            string filename = "exercise_b/" + 
                    methodName + to_string(ds[d]) +"d_" 
                    + to_string(Ns[N]) + "p_energy_alpha.txt";

            double alpha = alphaStart;
            for(int a=0; alpha > alphaStop+alphaStep; a++){

                int numberOfDimensions  = ds[d];
                int numberOfParticles   = Ns[N];
                int numberOfSteps       = (int) 1e6;
                double omega            = 1.0;          // Oscillator frequency.
                double stepLength       = 0.5;          // Metropolis step length.
                double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
                double timeStep         = 1e-2;

                System* system = new System();
                system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
                system->setWaveFunction             (new SimpleGaussian(system, alpha));
                system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
                system->setEquilibrationFraction    (equilibration);
                system->setStepLength               (stepLength);
                system->setAnalytical               (analyticOrNot);
                system->setFileName                 (filename);
                system->runMetropolisSteps          (numberOfSteps, firstCriteria);
                system->getSampler()->printOutputToEnergyAlphaFile()

                alpha -= alphaStep;
                firstCriteria = 1;
            }
        }
    }

*/

    // -------------------------------------------------------------------------------------

    // comparing numerical vs analytical
/*
    ofstream myfile;
    
    int firstCriteria = 0;                  // Criteria to print header in outputfile
    bool analyticOrNot = true;
    double beta = 1;
    double alpha            = 0.5;
    int numberOfDimensions  = 2;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
    double timeStep         = 1e-2;

    string methodName = setMethodName(analyticOrNot);

    string filename = methodName;
    
    myfile.open (filename, ios::out | ios::trunc);
    myfile.close(); 

    for (int n = 0; n < 10; n++){

        clock_t start, end;
        // Recording the starting clock tick.
        start = clock();

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (analyticOrNot);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria);
                
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    
        myfile.open (filename, ios::out | ios::app);
        
        myfile << time_taken << "\n";
        myfile.close();
    }
    */

// -------------------------------------------------------------------------------------

/* exercise c 

Implementing importance sampling and comparing with former results. 
Will the number of accepted steps increase?
*/

/*
    bool analyticOrNot = true;

    string methodName = setMethodName(analyticOrNot);

    std::vector <int> Ns = {1, 10, 100, 500};
    std::vector <int> ds = {1, 2, 3};

    std::vector <double> dls = {1.0, 0.5, 0.1, 0.01};

    ofstream thisfile;

    thisfile.open (methodName + "_importance", ios::out | ios::trunc);
    thisfile << " Step lenght: \t # of particles: \t Dimension: \t acceptance [%]: \t energy: \n";
    thisfile.close();

    for (int dl = 0; dl < dls.size(); dl++ ){

        for (int d = 0; d < ds.size(); d++){

            for (int N = 0; N < Ns.size(); N++ ){

                string filename = "exercise_b/importance_" + 
                        methodName + to_string(ds[d]) +"d_" 
                        + to_string(Ns[N]) + "p_energy_alpha.txt";


                int numberOfDimensions  = ds[d];
                int numberOfParticles   = Ns[N];
                int numberOfSteps       = (int) 1e6;
                double stepLength       = dls[dl];          // Metropolis step length.

                System* system = new System();
                system->setHamiltonian                (new HarmonicOscillator(system, omega, beta));
                system->setWaveFunction               (new SimpleGaussian(system, alpha));
                system->setInitialState               (new RandomUniform(system, numberOfDimensions, numberOfParticles));
                system->setEquilibrationFraction      (equilibration);
                system->setStepLength                 (stepLength);
                system->setAnalytical                 (analyticOrNot);
                system->setFileName                   (filename);
                system->runMetropolisStepsImportance  (numberOfSteps, firstCriteria);
                
                thisfile.open ("Output/exercis_c/" + methodName + "_importance.txt", ios::out | ios::app);
    
                thisfile << dls[dl] << "\t" << 
                            Ns[N] << "\t" << 
                            ds[d] << "\t" << 
                            system->getSampler()->getAcceptance() << "\t" <<
                            system->getSampler()->getEnergy() << "\n"; 
                thisfile.close();

            }
        }
    }
*/

/* exercise d

Need parallellization and also writing every energy to file

*/

    bool analyticOrNot = true;

    string methodName = setMethodName(analyticOrNot);

    numberOfDimensions    = 3;
    numberOfSteps         = (int) 1e6;
    stepLength            = 0.1;
    std::vector <int> Ns  = {1, 10, 100, 500};

    for (int N = 0; N < Ns.size(); N++ ){

        string filename = "exercise_d/importance_" + 
                methodName + to_string(numberOfDimensions) +"d_" 
                + to_string(Ns[N]) + "p_energy_alpha.txt";

        int numberOfParticles   = Ns[N];

        System* system = new System();
        system->setHamiltonian                (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction               (new SimpleGaussian(system, alpha));
        system->setInitialState               (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction      (equilibration);
        system->setStepLength                 (stepLength);
        system->setAnalytical                 (analyticOrNot);
        system->setFileName                   (filename);
        system->runMetropolisStepsImportance  (numberOfSteps, firstCriteria);

    }


}


void gradientDecentRun(string filename, 
                    double alphaGuess, 
                    double minimizationRate, 
                    double stopCriteria, 
                    bool analyticOrNot, 
                    bool importanceOrNot, 
                    bool interactionOrNot, 
                    double hardCoreDiameter, 
                    double beta){

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
        int numberOfParticles   = 2;
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double stepLength       = 0.5;          // Metropolis step length.
        double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
        double timeStep         = 1e-2;
        bool printToFileOrNot   = false;
        double energy           = 0;
    

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInteractionOrNot         (interactionOrNot, hardCoreDiameter); // Must come before initialize state
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (analyticOrNot);
        system->setImportanceSampling       (importanceOrNot, timeStep);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria);
    
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
                double hardCoreDiameter,
                double beta){

    int firstCriteria = 0;                  // Criteria to print header in outputfile

    double alpha = alphaStart;
    for(int a=0; alpha > alphaStop+alphaStep; a++){

        double stopCriteria     = 1e-9;         // Stopping criteria for energy vs exact energy.
        int numberOfDimensions  = 2;
        int numberOfParticles   = 2;
        int numberOfSteps       = (int) 1e6;
        double omega            = 1.0;          // Oscillator frequency.
        double stepLength       = 0.5;          // Metropolis step length.
        double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
        double timeStep         = 1e-2;
        bool printToFileOrNot   = true;

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setStepLength               (stepLength);
        system->setAnalytical               (analytic);
        system->setFileName                 (filename);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria);
            
        alpha -= alphaStep;
        firstCriteria = 1;
    }
 

}

string setMethodName(bool analyticOrNot){
    string methodName;

    if (analyticOrNot == true){
        methodName = "analytical_";
    }else{
        methodName = "numerical_";
    }
return methodName;
}