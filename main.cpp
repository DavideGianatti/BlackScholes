#include "BS.h"
using namespace std;

int main()
{
  double sigma = 0.2;   // market's volatility
  double r = 0.1;       // market's risk free interest rate
  double K = 50;        // strike price
  double T = 10;         // time to maturity

  int N_t = 16000;       // number of time steps
  int N_S = 200;        // number of grid points

  BS bs(sigma, r, K, T, N_t, N_S);
  bs.step(N_t);
  bs.save("bs.dat");
  cout << "Mean absolute error: " << bs.get_error_L1() << endl;

  BS bs_linearCC(sigma, r, K, T, N_t, N_S);
  bs_linearCC.step_linearCC(N_t);
  bs_linearCC.save("bs_linearCC.dat");
  cout << "Mean absolute error (with moving boundary conditions): " << bs_linearCC.get_error_L1() << endl;

  return 0;
}
