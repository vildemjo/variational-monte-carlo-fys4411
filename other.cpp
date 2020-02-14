#include <iostream>
#include <cmath>
#include <vector>
#include "random.h"
#include <iomanip>

// g++ -c random.cpp -o random.o
// g++ -o sample_new main.cpp random.o -std=c++0x

using namespace std;

double waveFunction(double alpha, double x){
    return exp(-alpha*x*x);
}

double waveFunctionMany(double alpha, vector<double> xs){
    double xSum;
    for(int n=0; n<xs.size(); n++){
        xSum += xs[n]*xs[n];
    }
    return exp(-alpha*xSum);
}

double localEnergyAnalytical(double alpha, double x){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;
    return (-hbar*hbar/(2*m))*(-2*alpha+4*x*x*alpha*alpha) + 0.5*m*omega*x*x;
}

double localEnergyAnalyticalMany(double alpha, vector <double> xs){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double numberOfParticles = xs.size();

    double xSum2 = 0.0;

    for(int n=0; n<numberOfParticles; n++){
        xSum2 += xs[n]*xs[n];
    }
    return (-hbar*hbar/(2*m))*(-2*alpha*numberOfParticles + 4*alpha*alpha*xSum2) + 0.5*m*omega*xSum2;
}

double localEnergyNumerical(double alpha, double x){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double energy;
    double dx = 0.001;
    double dpsidx2 = 0.0;
    
    dpsidx2 = (waveFunction(alpha,x+dx)+waveFunction(alpha,x-dx)-2*waveFunction(alpha,x))/(dx*dx);
    energy = (1.0/waveFunction(alpha, x))*(-hbar*hbar/(2*m)*dpsidx2) + 0.5*m*omega*x*x;

    return energy;
}

double localEnergyNumericalMany(double alpha, vector <double> xs){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    int length = xs.size();

    double energy = 0.0;
    double dx = 0.001;
    double dpsidx2 = 0.0;
    double xSum2 = 0.0;

    vector<double> x_less(length, 0.0);
    vector<double> x_more(length, 0.0);

    for (int n=0; n<length; n++){
        x_less[n] = xs[n]-dx;
        x_more[n] = xs[n]+dx;

        xSum2 += xs[n]*xs[n];
    }

    dpsidx2 = (waveFunctionMany(alpha,x_more)+waveFunctionMany(alpha,x_less)-2*waveFunctionMany(alpha,xs))/(dx*dx);
    energy += (1.0/waveFunctionMany(alpha, xs))*(-hbar*hbar/(2*m)*dpsidx2) + 0.5*m*omega*xSum2;
    
    return energy;
}

void MonteCarloSamplingOne(int maxVariations){
    Random random;

    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    int numberOfMCcycles = 1e6;
    double stepSize = 0.001;
    double positionOld = 0.0;
    double positionNew = 0.0;

    double waveFunctionOld = 0.0;
    double waveFunctionNew = 0.0;

    double deltaEnergy_ana = 0.0;
    double deltaEnergy_num = 0.0;

    double variance = 0.0;
    double error = 0.0;

    double alpha = 0.2;

    //double exactEnergy = k*k*hbar/(2*m);

    //cout << "Exact energy:" << exactEnergy <<  endl;
    cout << "alpha" << "\t" <<"energy (ana): " << "\t" << "energy (num): " << "\t" << "variance (ana): " << "\t" << "position: " << endl;

    // Evaluating the energy for different alpha (i.e. different wavefunctions)
    for (int i=0; i < maxVariations; i++){
        alpha += 0.05;
        
        double energy_ana = 0.0;
        double energy2_ana = 0.0;
        double energy_num = 0.0;
        double energy2_num = 0.0;

        // Pick a position:
        positionOld = stepSize*(random.nextDouble()-0.5);
        waveFunctionOld = waveFunction(alpha, positionOld);

        
        for (int j=0; j<numberOfMCcycles; j++){

            // Suggest a move to a new position:
            positionNew = positionOld + stepSize*(random.nextDouble()-0.5);
            waveFunctionNew = waveFunction(alpha, positionNew);

            // Check if new position is accepted:
            if (random.nextDouble() <= (waveFunctionOld*waveFunctionOld)/(waveFunctionNew*waveFunctionNew)){
                // If accepted, make the move:
                positionOld = positionNew;
                waveFunctionOld = waveFunctionNew;
            }
            // Calculate the local energy (in position x with alpha) analytically and numerically:
            deltaEnergy_ana = localEnergyAnalytical(alpha, positionOld);
            deltaEnergy_num = localEnergyNumerical(alpha, positionOld);
            // Calculating the sum of all local energies:
            energy_ana += deltaEnergy_ana;
            energy2_ana += deltaEnergy_ana*deltaEnergy_ana;

            energy_num += deltaEnergy_num;
            energy2_num += deltaEnergy_num*deltaEnergy_num; 
        }
        // Evaluating the expectation value - the average of all local energies:
        energy_ana /= numberOfMCcycles;
        energy2_ana /= numberOfMCcycles;
        energy_num /= numberOfMCcycles;
        energy2_num /= numberOfMCcycles;

        // Calculating the variance:
        variance = energy2_ana - energy_ana*energy_ana;
        
        // Estimating an error:
        error = sqrt(variance/numberOfMCcycles);
        cout << alpha << "\t" << energy_ana << "\t" <<  energy_num << "\t" << variance << "\t" << positionOld << endl;

    }
}

void MonteCarloSampling(int maxVariations, int numberOfParticles, double stepSize){
    Random random;

    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    int numberOfMCcycles = 1e6;

    vector<double> positionOld(numberOfParticles, 0.0);
    vector<double> positionNew(numberOfParticles, 0.0);

    double waveFunctionOld = 0.0;
    double waveFunctionNew = 0.0;

    double deltaEnergy_ana = 0.0;
    double deltaEnergy_num = 0.0;

    double variance = 0.0;
    double error = 0.0;

    double alpha = 0.3;

    cout << "alpha" << "\t" <<"energy (ana): " << "\t" << "energy (num): " << "\t" << "variance (ana): " << "\t" << "error: " << endl;

    // Evaluating the energy for different alpha (i.e. different wavefunctions)
    for (int i=0; i < maxVariations; i++){
        alpha += 0.05;
        
        double energy_ana = 0.0;
        double energy2_ana = 0.0;
        double energy_num = 0.0;
        double energy2_num = 0.0;

        // Pick a position:
        for(int m=0; m<numberOfParticles;m++){
            positionOld[m] = stepSize*(random.nextDouble()-0.5);
        }

        waveFunctionOld = waveFunctionMany(alpha, positionOld);

        // Moving one particle at the time
        for(int mm=0; mm<numberOfParticles; mm++){
            
            for (int j=0; j<numberOfMCcycles; j++){

                // Suggest a move to a new position:
                positionNew[mm] = positionOld[mm] + stepSize*(random.nextDouble()-0.5);
            
                waveFunctionNew = waveFunctionMany(alpha, positionNew);

                // Check if new position is accepted:
                if (random.nextDouble() <= (waveFunctionNew*waveFunctionNew)/(waveFunctionOld*waveFunctionOld)){
                    // If accepted, make the move:
                    positionOld = positionNew;
                    waveFunctionOld = waveFunctionNew;
                }
                // Calculate the local energy (in position x with alpha) analytically and numerically:
                deltaEnergy_ana = localEnergyAnalyticalMany(alpha, positionOld);
                deltaEnergy_num = localEnergyNumericalMany(alpha, positionOld);
                // Calculating the sum of all local energies:
                energy_ana += deltaEnergy_ana;
                energy2_ana += deltaEnergy_ana*deltaEnergy_ana;

                energy_num += deltaEnergy_num;
                energy2_num += deltaEnergy_num*deltaEnergy_num; 
            }
        }
        // Evaluating the expectation value - the average of all local energies:
        energy_ana /= (numberOfMCcycles*numberOfParticles);
        energy2_ana /= (numberOfMCcycles*numberOfParticles);
        energy_num /= (numberOfMCcycles*numberOfParticles);
        energy2_num /= (numberOfMCcycles*numberOfParticles);

        // Calculating the variance:
        variance = energy2_ana - energy_ana*energy_ana;
        
        // Estimating an error:
        error = sqrt(variance/numberOfMCcycles);
        cout << alpha << "\t" << energy_ana << "\t" << energy_num << "\t"  << variance << "\t" << error << endl;

    }
}


int main() {
    Random random;

    int maxVariantions = 10;
    int numberOfParticles = 1;
    double stepSize = 1;

    std::cout << std::setprecision(6) << std::fixed;

    MonteCarloSampling(maxVariantions, numberOfParticles, stepSize);
    // MonteCarloSamplingOne(maxVariantions);



    return 0;
}