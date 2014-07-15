#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "pimc.h"
#include "path.h"

using namespace std;

class path;

pimc::pimc(path p)
{
  observableSkip=50;
  equilSkip=1000;
  cnumAccept=0;// count accepted center of mass update
  snumAccept=0;// count accepted staging update
  vector<double>  EnergyTrace();
  
  //cout<<observableSkip<<endl;
  //cout<<equilSkip<<endl;
}

bool pimc::CenterOfMassMove(path p, int ptcl)
{
  srand((unsigned)time(NULL));

  int NumTime=p.get_NumTime();
  double delta = 0.5;
  double shift = delta*(-1.0 + 2.0*((double) rand() / (RAND_MAX)));
  double oldbeads[NumTime];

  // Store the positions on the worldline
  for (int slice=0; slice<NumTime; slice++)
  {
    oldbeads[slice] = p.beads[slice][ptcl];
  }

  // Calculate the potential action
  double oldAction = 0.0;
  for (int slice=0; slice<NumTime; slice++)
    oldAction += p.PotentialAction(slice);

  // Displace the worldline
  for (int slice=0; slice<NumTime; slice++)
    p.beads[slice][ptcl] = oldbeads[slice] + shift;

  // Compute the new action
  double newAction = 0.0;
  for (int slice=0; slice<NumTime; slice++)
    newAction += p.PotentialAction(slice);

  // Accept the move, or reject and restore the bead positions
  if( ((double) rand() / (RAND_MAX)) < exp(-(newAction - oldAction)) )
    return true;
  else
  {
    for (int slice=0; slice<NumTime; slice++)
    {
      p.beads[slice][ptcl] = oldbeads[slice];
    }
    return false;
  }
}

bool pimc::StagingMove(path p, int ptcl)
{
  // the length of the stage
  int m = 16;
  int NumTime=p.get_NumTime();
  double tau=p.get_tau();
  double lam=p.get_lam();
  double oldbeads[m-1];

  // Choose the start and end of the stage
  int alpha_start = ((int) NumTime*((double) rand() / (RAND_MAX)) ); //np.random.randint(0,Path.numTimeSlices)
  int alpha_end = (alpha_start + m) % NumTime;

  // Record the positions of the beads to be updated and store the action
  for(int a=0; a<m-1; a++)
    oldbeads[a] = 0;

  double oldAction = 0.0;
  for (int a=1; a<m; a++)
  {
      int slice = (alpha_start + a) % NumTime;
      oldbeads[a-1] = p.beads[slice][ptcl];
      oldAction += p.PotentialAction(slice);
  }

  // Generate new positions and accumulate the new action
  double newAction = 0.0;
  for (int a=1; a<m; a++)
  {
    int slice = (alpha_start + a) % NumTime;
    int slicem1 = (slice - 1) % NumTime;
    double tau1 = (m-a)*tau;
    double avex = (tau1*p.beads[slicem1][ptcl] +
            tau*p.beads[alpha_end][ptcl]) / (tau + tau1);
    double sigma2 = 2.0*lam / (1.0 / tau + 1.0 / tau1);
    p.beads[slice][ptcl] = avex + sqrt(sigma2)*((double) rand() / (RAND_MAX));
    newAction += p.PotentialAction(slice);
  }   

  // Perform the Metropolis step, if we rejct, revert the worldline
    if( ((double) rand() / (RAND_MAX)) < exp(-(newAction - oldAction)) )
        return true;
    else
        for(int a=1; a<m; a++)
        {
            int slice = (alpha_start + a) % NumTime;
            p.beads[slice][ptcl] = oldbeads[a-1];
            return false;
        }
}

vector<double> pimc::pimc_run(int NumSteps, path p)
{
  int NumParticle=p.get_NumParticle();
 
  for (int steps=0; steps<NumSteps; steps++)
  {
      // for each particle try a center-of-mass random move
      for (int i=0; i<NumParticle; i++) 
      {
        int ptcl=( (int) NumParticle*((double) rand() / (RAND_MAX)) );         
        if( CenterOfMassMove(p,ptcl) )
            cnumAccept += 1;
      }

      // for each particle try a staging move
      for (int i=0; i<NumParticle; i++)
      {
        int ptcl=( (int) NumParticle*((double) rand() / (RAND_MAX)) );
        if ( StagingMove(p,ptcl))
              snumAccept += 1;
      }

      // measure the energy
      if ( (steps % observableSkip == 0) && (steps > equilSkip) )
      {
          double E = p.Energy();
          EnergyTrace.push_back(E);
      }
  }

  cout<< "Acceptance Ratios:" <<endl;
  cout<< "Center of Mass:   " << ((1.0*cnumAccept)/(NumSteps*NumParticle)) << endl;
  cout<< "Staging:          " << ((1.0*snumAccept)/(NumSteps*NumParticle)) << endl;
  return EnergyTrace;

}
