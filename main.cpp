#include <iostream>
#include <fstream>
#include "system.h"
#include "particle.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/simplegaussianinteraction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/ellipticalharmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/gaussiandistribution.h"
#include "Math/random.h"
#include <string>
#include <omp.h>
#include <cmath>

using namespace std;

void alphaListRun(string filename, int numberOfP);
string setMethodName(bool analyticOrNot);

int main() {

// omp_set_num_threads(4);
// int thread_num = omp_get_max_threads ();

/* The standard set-up */

    bool analyticOrNot       = true;
    bool importanceOrNot     = true;
    bool allEnergiesOrNot    = true;
    int equilibration        = 1e5;          // Number of the total steps used for equilibration
    int numberOfDimensions   = 3;
    int numberOfParticles    = 10;
    int numberOfSteps        = (int) pow(2.0,20.0);
    double omega             = 1.0;          // Oscillator frequency.
    double stepLength        = 0.05;          // Metropolis step length.
    int firstCriteria        = 0;            // print header in file
    double alpha             = 0.483195;
    double inititalizingStep = stepLength;

    double hardCoreDiameter  = 0.0043;
    int numberOfBins = 800;
    double densityLength = 10.0;

    // elliptical or spherical trap (2.82843 or 1.0)
    double beta = 2.82843;    // omega_normal^2/omega_ho^2

    /* Set-up to run and save local energies for every step to file*/

/*
    allEnergiesOrNot    = true;

    clock_t start, end;
    // Recording the starting clock tick.
    start = clock();


    System* system = new System();
    system->setHamiltonian                (new HarmonicOscillator(system, omega));
    system->setWaveFunction               (new SimpleGaussian(system, alpha));
    system->setInitialState               (new GaussianDistribution(system, numberOfDimensions, 
                                                numberOfParticles, inititalizingStep));
    system->setEquilibration              (equilibration);
    system->setAnalytical                 (analyticOrNot);
    system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    system->setFileName                   ("Output/exercise_f/50p_3d_post_gradient_descent_normal");
    system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
                                                            allEnergiesOrNot, stepLength);

    cout << "number of steps: " << numberOfSteps << endl;

    cout << "energy: " << system->getSampler()->getEnergy()/((double) numberOfParticles*numberOfDimensions) << endl;
    cout << "should be: (0.504672)" << 0.5*alpha + 1/(8.0*alpha) << endl;

    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "CPU time: " << time_taken << " seconds" << endl;

*/

    /* Same but for the interacting case */   

/*
    importanceOrNot = false;
    stepLength = 0.5;
    inititalizingStep = stepLength;

    hardCoreDiameter  = 0.0043;
    numberOfBins = 800;
    densityLength = 10.0;

    // elliptical or spherical trap (2.82843 or 1.0)
    beta = 2.82843;    // omega_normal^2/omega_ho^2

    clock_t start, end;
    // Recording the starting clock tick.
    start = clock();


    System* system = new System();
    system->setHamiltonian                (new EllipticalHarmonicOscillator(system, omega, beta));
    system->setWaveFunction               (new SimpleGaussianInteraction(system, alpha, hardCoreDiameter, beta));
    system->setInitialState               (new RandomUniform(system, numberOfDimensions, 
                                                numberOfParticles, inititalizingStep));
    system->setEquilibration              (equilibration);
    system->setAnalytical                 (analyticOrNot);
    system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
    system->setFileName                   ("Output/exercise_f/50p_3d_post_gradient_descent_interacting");
    system->runMetropolisSteps            (numberOfSteps, firstCriteria, importanceOrNot, 
                                                            allEnergiesOrNot, stepLength);

    cout << "energy: " << system->getSampler()->getEnergy()/((double) numberOfParticles*numberOfDimensions) << endl;

    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "CPU time: " << time_taken << " seconds" << endl;

/*


/* Exercise b */

/* comparing with exact answear */

/*
    double alphaStart = 0.7;
    double alphaStop = 0.01;
    double alphaStep = 0.02;

    string methodName = setMethodName(analyticOrNot);

    std::vector <int> Ns = {1, 10, 50, 100};//{100, 500};
    std::vector <int> ds = {1};//{1, 2, 3};

    for (int d = 0; d < ds.size(); d++){

        for (int N = 0; N < Ns.size(); N++ ){

            ofstream energyfile;
            
            string filename = "exercise_b/" + 
                    methodName + to_string(ds[d]) +"d_" 
                    + to_string(Ns[N]) + "p_energy_alpha.txt";

            energyfile.open ("Output/" + filename, ios::out | ios::trunc);
            energyfile << "Alpha: \t Energy: \n";
            energyfile.close();

            alpha = alphaStart;
            for(int a=0; alpha > alphaStop+alphaStep; a++){

                double energy = 0;

                numberOfDimensions  = ds[d];
                numberOfParticles   = Ns[N];
                numberOfSteps       = (int) std::pow(2.0,19);
                allEnergiesOrNot    = false;
                
                System* system = new System();
                system->setHamiltonian              (new HarmonicOscillator(system, omega));
                system->setWaveFunction             (new SimpleGaussian(system, alpha));
                system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, inititalizingStep));
                system->setEquilibration            (equilibration);
                system->setAnalytical               (analyticOrNot);
                system->setFileName                 (filename);
                system->runMetropolisSteps          (numberOfSteps, firstCriteria, importanceOrNot, allEnergiesOrNot, stepLength);
            
                energy = system->getSampler()->getEnergy();

                energyfile.open ("Output/"+ filename, ios::out | ios::app);
                energyfile << alpha << "\t" << energy << "\n";
                energyfile.close();

                alpha -= alphaStep;
                firstCriteria = 1;
            }
        
        }
    }
*/

/* Printing all local energy measurements for some alphas */

/*
    double alphaStart = 0.65;
    double alphaStop = 0.35;
    double alphaStep = 0.05;
    
    allEnergiesOrNot = true;
    numberOfSteps = pow(2,20);

    string methodName = setMethodName(analyticOrNot);

    std::vector <int> Ns = {1, 10};//, };
    std::vector <int> ds = {3};

    for (int d = 0; d < ds.size(); d++){

        for (int N = 0; N < Ns.size(); N++ ){

            alpha = alphaStart;
            for(int a=0; alpha >= alphaStop; a++){

                double alphaPrint = alpha*100.0;
                int alphaPrintable = ceil(alphaPrint);    
                string filename = "Output/exercise_b/allEnergies/" + 
                    methodName + to_string(ds[d]) +"d_" 
                    + to_string(Ns[N])+"p"+ "_alpha_" + to_string(alphaPrintable) + "_MC_20";

                numberOfDimensions  = ds[d];
                numberOfParticles   = Ns[N];
                
                System* system = new System();
                system->setHamiltonian              (new HarmonicOscillator(system, omega));
                system->setWaveFunction             (new SimpleGaussian(system, alpha));
                system->setInitialState             (new RandomUniform(system, numberOfDimensions, 
                                                            numberOfParticles, inititalizingStep));
                system->setEquilibration            (equilibration);
                system->setAnalytical               (analyticOrNot);
                system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
                system->setFileName                 (filename);
                system->runMetropolisSteps          (numberOfSteps, firstCriteria, importanceOrNot, 
                                                                    allEnergiesOrNot, stepLength);
                
                alpha -= alphaStep;
                firstCriteria = 1;
            }
        }
    }

*/


/* comparing numerical vs analytical
(run once for "analyticalOrNot = true" and once "= false") */

/*
    ofstream myfile;
    
    analyticOrNot = true;
    numberOfDimensions  = 3;
    numberOfSteps       = (int) std::pow(2.0,19);
    allEnergiesOrNot = false;
    

    for (int m = 1; m <= 1; m +=2){ 

        numberOfParticles   = m;

        string methodName = setMethodName(analyticOrNot);

        string filename = "Output/exercise_b/" + methodName + to_string(m) + "p_3d_cpu_time.txt";
        
        myfile.open(filename, ios::out | ios::trunc);
        myfile.close(); 

        for (int n = 0; n < 10; n++){

            clock_t start, end;
            // Recording the starting clock tick.
            start = clock();

            System* system = new System();
            system->setHamiltonian              (new HarmonicOscillator(system, omega));
            system->setWaveFunction             (new SimpleGaussian(system, alpha));
            system->setInitialState             (new RandomUniform(system, numberOfDimensions, 
                                                        numberOfParticles, inititalizingStep));
            system->setEquilibrationFraction    (equilibration);
            system->setStepLength               (stepLength);
            system->setAnalytical               (analyticOrNot);
            system->runMetropolisSteps          (numberOfSteps, firstCriteria, importanceOrNot, 
                                                                    allEnergiesOrNot, stepLength);
                    

            end = clock();

            double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
        
            myfile.open (filename, ios::out | ios::app);
            
            myfile << time_taken << "\n";
            myfile.close();
        }
    }
*/
// -------------------------------------------------------------------------------------

/* exercise c 

Implementing importance sampling and comparing with former results. 
*/
/*

    analyticOrNot    = true;
    allEnergiesOrNot = false;
    importanceOrNot  = false;

    numberOfDimensions  = 3;
    numberOfParticles   = 3;
    alpha               = 0.45;
    numberOfSteps       = (int) pow(2.0, 20.0);

    string methodName = setMethodName(analyticOrNot);
    string samplingType;

    std::vector <double> dls = {1.0, 0.5, 0.1, 0.05, 0.01, 0.005};

    ofstream thisfile;

    if (importanceOrNot == true){
        samplingType = "importance";
    }else{
        samplingType = "brute_force";
    }

    string file_name = "Output/exercise_c/" + methodName + 
                        to_string(numberOfParticles) + "p_" +
                        to_string(numberOfDimensions) + 
                        "d_" + samplingType + ".txt";

    thisfile.open(file_name, ios::out | ios::trunc);
    thisfile << " Step lenght: \t acceptance [%]: \t energy: \t CPU time: \n";
    thisfile.close();

    for (int dl = 0; dl < dls.size(); dl++ ){

        clock_t start, end;
        // Recording the starting clock tick.
        start = clock();

        double timeStep     = dls[dl];          // Metropolis step length.
        inititalizingStep = timeStep;

        System* system = new System();
        system->setHamiltonian                (new HarmonicOscillator(system, omega));
        system->setWaveFunction               (new SimpleGaussian(system, alpha));
        if (importanceOrNot == true){
            system->setInitialState           (new GaussianDistribution(system, numberOfDimensions, 
                                                        numberOfParticles, inititalizingStep));
        }else{
            system->setInitialState           (new RandomUniform(system, numberOfDimensions, 
                                                        numberOfParticles, inititalizingStep));
        }
        system->setEquilibration              (equilibration);
        system->setAnalytical                 (analyticOrNot);
        system->runMetropolisSteps            (numberOfSteps, firstCriteria, 
                                                importanceOrNot, allEnergiesOrNot, timeStep);
        
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 

        thisfile.open(file_name, ios::out | ios::app);

        thisfile << dls[dl] << "\t" << 
                    system->getSampler()->getAcceptance() << "\t" <<
                    system->getSampler()->getEnergy() << "\t"<<
                    time_taken << "\n"; 
        thisfile.close();

    }
*/
/* Comparing for different amounts of MC cycles */

/*
    analyticOrNot    = true;
    allEnergiesOrNot = false;
    importanceOrNot  = false;

    numberOfDimensions  = 3;
    numberOfParticles   = 3;
    alpha               = 0.45;

    string methodName = setMethodName(analyticOrNot);
    string samplingType;

    std::vector <double> MC = {10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0};

    ofstream thisfile;

    double timeStep;


    if (importanceOrNot == true){
        samplingType = "importance";
        timeStep = 0.005;
    }else{
        samplingType = "brute_force";
        timeStep = 0.5;
    }

    string file_name = "Output/exercise_c/" + methodName + 
                        to_string(numberOfParticles) + "p_" +
                        to_string(numberOfDimensions) + 
                        "d_" + samplingType + "MC_cycles.txt";

    thisfile.open(file_name, ios::out | ios::trunc);
    thisfile << " Step lenght: \t acceptance [%]: \t energy: \t CPU time: \n";
    thisfile.close();

    for (int mc = 0; mc < MC.size(); mc++ ){

        numberOfSteps       = (int) pow(2.0, MC[mc]);

        clock_t start, end;
        // Recording the starting clock tick.
        start = clock();

                  // Metropolis step length.
        inititalizingStep = timeStep;

        System* system = new System();
        system->setHamiltonian                (new HarmonicOscillator(system, omega));
        system->setWaveFunction               (new SimpleGaussian(system, alpha));
        if (importanceOrNot == true){
            system->setInitialState           (new GaussianDistribution(system, numberOfDimensions, 
                                                        numberOfParticles, inititalizingStep));
        }else{
            system->setInitialState           (new RandomUniform(system, numberOfDimensions, 
                                                        numberOfParticles, inititalizingStep));
        }
        system->setEquilibration              (equilibration);
        system->setAnalytical                 (analyticOrNot);
        system->runMetropolisSteps            (numberOfSteps, firstCriteria, 
                                                importanceOrNot, allEnergiesOrNot, timeStep);
        
        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 

        thisfile.open(file_name, ios::out | ios::app);

        thisfile << MC[mc] << "\t" << 
                    system->getSampler()->getAcceptance() << "\t" <<
                    system->getSampler()->getEnergy() << "\t"<<
                    time_taken << "\n"; 
        thisfile.close();

    }
*/

/* Comparing with exercise b */

/*

    double alphaStart = 0.65;
    double alphaStop = 0.35;
    double alphaStep = 0.05;
    importanceOrNot = true;
    stepLength = 0.005;
    inititalizingStep = stepLength;

    allEnergiesOrNot = true;
    numberOfSteps = pow(2,20);

    string methodName = setMethodName(analyticOrNot);

    std::vector <int> Ns = {1, 10};//, };
    std::vector <int> ds = {3};

    for (int d = 0; d < ds.size(); d++){

        for (int N = 0; N < Ns.size(); N++ ){

            alpha = alphaStart;
            for(int a=0; alpha >= alphaStop; a++){

                double alphaPrint = alpha*100.0;
                int alphaPrintable = ceil(alphaPrint);    
                string filename = "Output/exercise_c/allEnergies/" + 
                    methodName + to_string(ds[d]) +"d_" 
                    + to_string(Ns[N])+"p"+ "_alpha_" + to_string(alphaPrintable) + "_MC_20";

                numberOfDimensions  = ds[d];
                numberOfParticles   = Ns[N];
                
                System* system = new System();
                system->setHamiltonian              (new HarmonicOscillator(system, omega));
                system->setWaveFunction             (new SimpleGaussian(system, alpha));
                system->setInitialState             (new GaussianDistribution(system, numberOfDimensions, 
                                                            numberOfParticles, inititalizingStep));
                system->setEquilibration            (equilibration);
                system->setAnalytical               (analyticOrNot);
                system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
                system->setFileName                 (filename);
                system->runMetropolisSteps          (numberOfSteps, firstCriteria, importanceOrNot, 
                                                                    allEnergiesOrNot, stepLength);
                
                alpha -= alphaStep;
                firstCriteria = 1;
            }
        }
    }
*/

// ----------------------------------------------------------------

/*  Exercise e 

    Implemting interaction
*/


/* comparing with b and c */

/*
    double alphaStart = 0.65;
    double alphaStop = 0.35;
    double alphaStep = 0.05;
    importanceOrNot = false;
    stepLength = 0.5;
    inititalizingStep = stepLength;

    hardCoreDiameter  = 0.0043;
    numberOfBins = 800;
    densityLength = 10.0;

    // elliptical or spherical trap (2.82843 or 1.0)
    beta = 2.82843;    // omega_normal^2/omega_ho^2

    allEnergiesOrNot = true;
    numberOfSteps = pow(2,20);

    string methodName = setMethodName(analyticOrNot);

    std::vector <int> Ns = {50};//, };
    std::vector <int> ds = {3};

    for (int d = 0; d < ds.size(); d++){

        for (int N = 0; N < Ns.size(); N++ ){

            alpha = alphaStart;
            for(int a=0; alpha > alphaStop; a++){

                double alphaPrint = alpha*100.0;
                int alphaPrintable = ceil(alphaPrint);    
                string filename = "Output/exercise_e/allEnergies/" + 
                    methodName + to_string(ds[d]) +"d_" 
                    + to_string(Ns[N])+"p"+ "_alpha_" + to_string(alphaPrintable) + "_MC_20";

                numberOfDimensions  = ds[d];
                numberOfParticles   = Ns[N];
                
                System* system = new System();
                system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega, beta));
                system->setWaveFunction             (new SimpleGaussianInteraction(system, alpha, hardCoreDiameter, beta));
                system->setInitialState             (new RandomUniform(system, numberOfDimensions, 
                                                            numberOfParticles, inititalizingStep));
                system->setEquilibration            (equilibration);
                system->setAnalytical               (analyticOrNot);
                system->getWaveFunction               ()->setOneBodyDensityBins(numberOfBins, densityLength);
                system->setFileName                 (filename);
                system->runMetropolisSteps          (numberOfSteps, firstCriteria, importanceOrNot, 
                                                                    allEnergiesOrNot, stepLength);
                
                alpha -= alphaStep;
                firstCriteria = 1;
            }
        }
    }
*/

// ---------------------------------------------------------------

/* exercise f

Optimization with use of gradient descent
*/

/*
    ofstream file;
    
    double energyChange = 1.0;
    double stopCriteria = 1e-7;
    double energyNew = 0.0;
    double energyDerivative = 1.0;
    double alphaNew = 0;
    double minimizationRate = 0.1;
    allEnergiesOrNot = false;
    importanceOrNot = true;
    alpha = 0.51;
    stepLength = 0.005;
    inititalizingStep = stepLength;

    double energy       = 0;

    numberOfDimensions  = 3;
    numberOfParticles   = 10;
    numberOfSteps       = (int) std::pow(2,19.0);

    string file_name = "Output/exercise_f/m_1_10p_3d_alpha_51_gradient_descent_normal.txt";

    file.open (file_name, ios::out | ios::trunc);
    file << "Alpha: \t Energy: \t Derivative: \n";
    file.close();


    for (int k=0;  energyChange > stopCriteria || k == 60; k++){

        clock_t start, end;
    // Recording the starting clock tick.
        start = clock();
    

        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new GaussianDistribution(system, numberOfDimensions, 
                                                    numberOfParticles, inititalizingStep));
        system->setEquilibration            (equilibration);
        system->setAnalytical               (analyticOrNot);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria, 
                                            importanceOrNot, allEnergiesOrNot, stepLength);

        firstCriteria = 1;
        
        energyNew = system->getSampler()->getEnergy();
        energyDerivative = system->getSampler()->getDerivative();
        alphaNew = alpha - minimizationRate*energyDerivative/numberOfParticles;

        
        file.open (file_name, ios::out | ios::app);
        file << alpha << "\t" << energy << "\t" << energyDerivative << "\n";
        file.close();

        energyChange = std::abs(energyNew - energy);
        alpha = alphaNew;
        energy = energyNew;

        end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
        cout << "CPU time: " << time_taken << " seconds" << endl;
    } 
*/

    /* A little more complex gradient descent */

/*
    ofstream file;
    
    double energyChange = 1.0;
    double stopCriteria = 1e-5;
    double energyNew = 0.0;
    double energyDerivative = 0;
    double alphaNew = 0;
    stepLength = 0.005;
    inititalizingStep = stepLength;
    importanceOrNot = true;
    allEnergiesOrNot = false;

    double minimizationRate = 0.1;
    double previousDerivativeWeight = 0.03;
    double energyDerivativePrevious = 0;
    double energy = 0;
    alpha = 0.65;

    numberOfDimensions  = 3;
    numberOfParticles   = 10;
    numberOfSteps       = (int) std::pow(2,19.0);

    string nameOfFile = "exercise_f/m_1_w_03_" + to_string(numberOfDimensions) + "d_"
                        + to_string(numberOfParticles) + "p_alpha_65_gradient_descent_extended";

    file.open ("Output/" + nameOfFile + ".txt", ios::out | ios::trunc);
    file << "Alpha: \t Energy: \t Derivative: \n";
    file.close();


    for (int k=0;  energyChange > stopCriteria; k++){


        System* system = new System();
        system->setHamiltonian              (new HarmonicOscillator(system, omega));
        system->setWaveFunction             (new SimpleGaussian(system, alpha));
        system->setInitialState             (new GaussianDistribution(system, numberOfDimensions,
                                                     numberOfParticles, inititalizingStep));
        system->setEquilibration            (equilibration);
        system->setAnalytical               (analyticOrNot);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria, importanceOrNot,
                                                         allEnergiesOrNot, stepLength);

        firstCriteria = 1;
        
        energyNew = system->getSampler()->getEnergy();
        energyDerivativePrevious = energyDerivative;
        energyDerivative = system->getSampler()->getDerivative();
        alphaNew = alpha - (minimizationRate-previousDerivativeWeight)*energyDerivative/numberOfParticles - previousDerivativeWeight*energyDerivativePrevious/numberOfParticles;
        
        file.open ("Output/"+ nameOfFile + ".txt", ios::out | ios::app);
        file << alpha << "\t" << energy << "\t" << energyDerivative << "\n";
        file.close();

        energyChange = std::abs(energyNew - energy);
        alpha = alphaNew;
        energy = energyNew;
    }    

*/

/* finding alpha for interacting case with gradient descent method */

/*
    ofstream file;
    
    double energyChange = 1.0;
    double stopCriteria = 1e-7;
    double energyNew = 0.0;
    double energyDerivative = 1.0;
    double alphaNew = 0;
    double minimizationRate = 0.1;
    allEnergiesOrNot = false;
    importanceOrNot = false;
    alpha = 0.5;
    stepLength = 0.5;
    inititalizingStep = stepLength;

    hardCoreDiameter  = 0.0043;
    numberOfBins = 800;
    densityLength = 10.0;

    // elliptical or spherical trap (2.82843 or 1.0)
    beta = 2.82843;    // omega_normal^2/omega_ho^2

    double energy       = 0;

    numberOfDimensions  = 3;
    numberOfParticles   = 10;
    numberOfSteps       = (int) std::pow(2,19.0);

    string file_name = "Output/exercise_f/m_1_10p_3d_interaction_no_jastrow.txt";

    file.open (file_name, ios::out | ios::trunc);
    file << "Alpha: \t Energy: \t Derivative: \n";
    file.close();


    for (int k=0;  energyChange > stopCriteria; k++){
    

        System* system = new System();
        system->setHamiltonian              (new EllipticalHarmonicOscillator(system, omega, beta));
        system->setWaveFunction             (new SimpleGaussianInteraction(system, alpha, hardCoreDiameter, beta));
        system->setInitialState             (new RandomUniform(system, numberOfDimensions, 
                                                    numberOfParticles, inititalizingStep));
        system->setEquilibration            (equilibration);
        system->setAnalytical               (analyticOrNot);
        system->runMetropolisSteps          (numberOfSteps, firstCriteria, 
                                            importanceOrNot, allEnergiesOrNot, stepLength);

        firstCriteria = 1;
        
        energyNew = system->getSampler()->getEnergy();
        energyDerivative = system->getSampler()->getDerivative();
        alphaNew = alpha - minimizationRate*energyDerivative/numberOfParticles;

        
        file.open (file_name, ios::out | ios::app);
        file << alpha << "\t" << energy << "\t" << energyDerivative << "\n";
        file.close();

        energyChange = std::abs(energyNew - energy);
        alpha = alphaNew;
        energy = energyNew;
    } 
*/


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