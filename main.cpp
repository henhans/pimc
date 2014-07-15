/* Pathe integral Monte Carlo code based on Dr. Adrian Del Maesstrc's lecture notes. I try
   to translate it to c++.                                Author:Tsung-Han Lee 7/13/2014 */

#include "path.h"
#include "pimc.h"
//#include "routines.cpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <numeric>

using namespace std;

int main(int argc, char* argv[])
{
  double T = 1.25;  // temperature in Kelvin  
  double lam = 0.5;    // \hbar^2/2m k_B

  int  numParticle = 1; // number of particles
  int  numTime = 20; //number of time slices
  int  numMCSteps = 100000;//00000; 

  double  tau = 1.0/(T*numTime);

  vector<double> Energy;

  //print informations
  cout<< "Simulation Parameters:" << endl;
  cout<< "N      = " << numParticle << endl;
  cout<< "tau    = " << tau << endl;
  cout<< "lambda = " << lam << endl;
  cout<< "T      = " << T   << endl;

  //srand((unsigned)time(NULL));
  srand(1234);

  path p(numParticle, numTime, tau, lam);

  //p.PotentialAction(19);
  pimc mc(p);

  Energy= mc.pimc_run( numMCSteps, p);

  double sum = accumulate(Energy.begin(), Energy.end(), 0.0);
  double mean = sum / Energy.size();

  double sq_sum = inner_product(Energy.begin(), Energy.end(), Energy.begin(), 0.0);
  double stdev = sqrt(sq_sum / Energy.size() - mean * mean)/sqrt(Energy.size());

  cout<< "Energy = "<< mean << " +/- "<< stdev  <<endl;

  for( vector<double>::const_iterator i=Energy.begin(); i!=Energy.end(); ++i)
     cout << *i<< " ";
  cout<<endl;

  return 0;
}
