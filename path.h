using namespace std;

class path
{
  public:
  path( int NP, int NT, double t, double l);

  double Vext(int N, double* x);
  double PotentialAction(int slice);
  double KineticEnergy();
  double PotentialEnergy();
  double Energy();
  int get_NumTime();
  int get_NumParticle();
  double get_tau();
  double get_lam();

  double** beads;

  private:
  int NumParticle;
  int NumTime;
  //const int N;
  double tau;
  double lam;
};

