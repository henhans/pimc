#include <vector>

using namespace std;

class path;

class pimc
{
  public:
  pimc(path p);

  vector<double> pimc_run(int NumSteps, path p);
  bool CenterOfMassMove(path p, int ptcl);
  bool StagingMove(path p, int ptcl);

  private:
  int observableSkip;
  int equilSkip;
  int cnumAccept;// count accepted center of mass update
  int snumAccept;// count accepted staging update
  vector<double>  EnergyTrace;
};
