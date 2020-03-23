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
#include <cmath>

using namespace std;

void gradientDecentRun(string filename, double alphaStart, double minimizationRate, double stopCriteria, bool analytic, bool importance, bool interactionOrNot, double hardCoreDiameter, double beta);
void alphaListRun(string filename, int numberOfP);
string setMethodName(bool analyticOrNot);

int main() {

// omp_set_num_threads(4);
// int thread_num = omp_get_max_threads ();

/* The standard set-up */

    bool analyticOrNot      = true;
    double hardCoreDiameter = 0.0043;
    int numberOfDimensions  = 3;
    int numberOfParticles   = 3;
    int numberOfSteps       = (int) pow(2.0,21.0);
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 1.0;        // Metropolis step length.
    double equilibration    = 0.0;          // Fraction of the total steps used for equilibration
    double timeStep         = 0.0003;
    int firstCriteria       = 0;            // print header in file
    double alpha            = 0.5;

    int numberOfBins = 200;
    double densityLength = 4.0;

    // elliptical or spherical trap (2.82843 or 1.0)
    double beta = 1.0; //2.82843;    // omega_normal^2/omega_ho^2


    System* system = new System();
    system->setHamiltonian                (new HarmonicOscillator(system, omega, beta));
    system->setWaveFunction               (new SimpleGaussian(system, alpha));
    system->setInitialState               (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction      (equilibration);
    system->setAnalytical                 (analyticOrNot);
    system->getWaveFunction()->setOneBodyDensityBins(numberOfBins, densityLength);
    system->setFileName                   ("Output/testing_density_23_03_");
    system->runMetropolisSteps            (numberOfSteps, firstCriteria, stepLength);


    cout << "number of steps: " << numberOfSteps << endl;
    



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

                System* system = new System();
                system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
                system->setWaveFunction             (new SimpleGaussian(system, alpha));
                system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
                system->setEquilibrationFraction    (equilibration);
                system->setAnalytical               (analyticOrNot);
                system->setFileName                 (filename);
                system->runMetropolisSteps          (numberOfSteps, firstCriteria, stepLength);
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
    
    analyticOrNot = true;
    numberOfDimensions  = 3;
    numberOfParticles   = 10;
    numberOfSteps       = (int) 1e6;
    stepLength          = 0.5;          // Metropolis step length.
    equilibration       = 0.1;          // Fraction of the total steps used for equilibration
    
    string methodName = setMethodName(analyticOrNot);

    string filename = "Output/" + methodName + ".txt";
    
    // myfile.open (filename, ios::out | ios::trunc);
    // myfile.close(); 

    for (int n = 0; n < 1; n++){

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
        system->runMetropolisSteps          (numberOfSteps, firstCriteria, stepLength);
                

        end = clock();

        double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    
        // myfile.open (filename, ios::out | ios::app);
        
        // myfile << time_taken << "\n";
        // myfile.close();
        std::cout << time_taken << std::endl;
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
                double timeStep         = dls[dl];          // Metropolis step length.

                System* system = new System();
                system->setHamiltonian                (new HarmonicOscillator(system, omega, beta));
                system->setWaveFunction               (new SimpleGaussian(system, alpha));
                system->setInitialState               (new RandomUniform(system, numberOfDimensions, numberOfParticles));
                system->setEquilibrationFraction      (equilibration);
                system->setAnalytical                 (analyticOrNot);
                system->setFileName                   (filename);
                system->runMetropolisStepsImportance  (numberOfSteps, firstCriteria, timeStep);
                
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

/*
    analyticOrNot = true;

    string methodName = setMethodName(analyticOrNot);

    numberOfDimensions    = 3;
    numberOfSteps         = (int) 1e6;
    timeStep              = 0.1;
    std::vector <int> Ns  = {1, 10, 100, 500};

    for (int N = 0; N < Ns.size(); N++ ){

        string filename = "exercise_d/importance_" + 
                methodName + to_string(numberOfDimensions) +"d_" 
                + to_string(Ns[N]) + "p";

        numberOfParticles   = Ns[N];

        System* system = new System();
        system->setHamiltonian                (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction               (new SimpleGaussian(system, alpha));
        system->setInitialState               (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction      (equilibration);
        system->setAnalytical                 (analyticOrNot);
        system->setFileName                   (filename);
        system->runMetropolisStepsImportance  (numberOfSteps, firstCriteria, timeStep);

    }

*/

// string filename = "new_test_2p_3d";
// for (int n = 1; n <= 5; n++){
// string filename = "total_energy_"+ to_string(n) + "p_1d";
// alphaListRun(filename, n);
// }
}


// void gradientDecentRun(string filename, 
//                     double alphaGuess, 
//                     double minimizationRate, 
//                     double stopCriteria, 
//                     bool analyticOrNot, 
//                     bool importanceOrNot, 
//                     bool interactionOrNot, 
//                     double hardCoreDiameter, 
//                     double beta){

//     int firstCriteria               = 0;                  // Criteria to print header in outputfile

//     double alpha                    = alphaGuess;         // Variational parameter first guess
//     double m_stopCriteria           = stopCriteria;       // Stopping criteria for energy vs exact energy.
//     double m_minimizationRate       = minimizationRate;   // Model for double derivative of energy.
//     double energyDerivative         = 1.0;
//     double alphaChange              = 0.5;
//     double energyChange             = 1;
//     double alphaNew                 = alpha;
//     double energyNew                = 0.8; 
//     double gradientPrintToFileOrNot = false;

//     ofstream file;
//     if (gradientPrintToFileOrNot == true){
//         file.open ("Output/gradient_descent.txt", ios::out | ios::trunc);
//         file << "Alpha: \t Energy: \t Derivative: \n";
//         file.close();
//     }

//     for (int k=0;  energyChange > m_stopCriteria; k++){

//         int numberOfDimensions  = 1;
//         int numberOfParticles   = 2;
//         int numberOfSteps       = (int) 1e6;
//         double omega            = 1.0;          // Oscillator frequency.
//         double stepLength       = 0.5;          // Metropolis step length.
//         double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
//         double timeStep         = 1e-2;
//         bool printToFileOrNot   = false;
//         double energy           = 0;
    

//         System* system = new System();
//         system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
//         system->setWaveFunction             (new SimpleGaussian(system, alpha));
//         system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
//         system->setEquilibrationFraction    (equilibration);
//         system->setStepLength               (stepLength);
//         system->setAnalytical               (analyticOrNot);
//         system->runMetropolisSteps          (numberOfSteps, firstCriteria);
    
//         firstCriteria = 1;
        
//         energyNew = system->getSampler()->getEnergy();
//         energyDerivative = system->getSampler()->getDerivative();
//         alphaNew = alpha - m_minimizationRate*energyDerivative/numberOfParticles;
        
//         energyChange = std::abs(energyNew - energy);
//         alpha = alphaNew;
//         energy = energyNew;

//         if (gradientPrintToFileOrNot == true){
//             file.open ("Output/gradient_descent.txt", ios::out | ios::app);
//             file << alpha << "\t" << energy << "\t" << energyDerivative << "\n";
//             file.close();
//         }
        
//     }
    
// }




void alphaListRun(string filename, int numberOfP){

    bool analyticOrNot      = true;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Fraction of the total steps used for equilibration
    int firstCriteria       = 0;            // print header in file

    // elliptical or spherical trap (2.82843 or 1.0)
    double beta = 1.0;//2.82843;    // omega_normal^2/omega_ho^2

    double alphaStart = 0.6;
    double alphaStop = 0.01;
    double alphaStep = 0.02;

    int numberOfDimensions  = 1;
    int numberOfParticles   = numberOfP;

    double alpha = alphaStart;
    for(int a=0; alpha > alphaStop+alphaStep; a++){

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega, beta));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction    (equilibration);
        system->setAnalytical               (analyticOrNot);
        system->setFileName                 (filename);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria, stepLength);
        
        // cout << "alpha: " << alpha << endl;
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