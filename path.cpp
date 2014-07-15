#include "path.h"
#include "routines.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <cmath>


using namespace std;

path::path( int NP, int NT, double t, double l)
{
  NumParticle = NP;
  NumTime = NT;
  tau = t;
  lam = l;

  //N = NumParticle;
  beads = new double*[NumTime];
  for( int i =0; i<NumTime; i++)
  {
    beads[i]= new double[NumParticle];
  }

  //initialize the random position for each  particle at each time slice
  //srand((unsigned)time(NULL)); // initialize random seed

  for(int i=0; i<NumTime; i++)
  {
    for( int j=0; j<NumParticle; j++)
    {
      beads[i][j] = 0.5*(-1.0 + 2.0*((double) rand() / (RAND_MAX))); 
    }
  }
  
  for(int i=0; i<NumTime; i++)
  {
    for( int j=0; j<NumParticle; j++)
    {
      //cout<< "bead["<< i <<"]["<<j<<"]= " << beads[i][j] <<endl;
    }
  }
 
}

double path::Vext(int N, double* x)
{
  return HarmonicOscillator(N, x);
}

double path::PotentialAction(int slice)
{
  double pot=0;

  pot=Vext(NumParticle, &(beads[slice][0]));
  //cout<< "s=" <<s<<"NumParticle"<< NumParticle <<endl;
  //cout<< "pot=" <<pot <<endl;
  return tau*pot;
}

double path::KineticEnergy()
{
  double tot = 0.0;
  double norm = 1.0/(4.0*lam*tau*tau);
  for (int slice=0; slice< NumTime; slice++)
  {
    int slicep1 = (slice + 1) % NumTime;//periodic boundary condition
    for (int ptcl=0; ptcl< NumParticle; ptcl++)
    {
        double delR = beads[slicep1][ptcl] - beads[slice][ptcl];
        tot = tot - norm*delR*delR;
    }
  }

  double KE = 0.5*NumParticle/tau + tot/(NumTime);
  return KE;
}

double path::PotentialEnergy()
{
  double PE = 0.0;
  for (int slice=0; slice<NumTime; slice++)
  {
    PE = PE + Vext(NumParticle, &beads[slice][0]);
  }

  return PE/(NumTime);

}

double path::Energy()
{
  return PotentialEnergy()+KineticEnergy();
}

int path::get_NumTime()
{
  return NumTime;
}

int path::get_NumParticle()
{
  return NumParticle;
}

double path::get_tau()
{
  return tau;
}

double path::get_lam()
{
  return lam;
}

