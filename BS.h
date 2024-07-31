#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>
using namespace std;

#define boundary 4   // boundary conditions set to boundary*K
#define S_tilda 1000000 // point fixed at infinity
#define see_till 2   // observe solution up to see_till * K = see_till * delta_S * N_S / boundary

/*
Implementation of an explicit method for solving the Black-Scholes parabolic PDE for
European Call options using finite differences.
*/

class BS
{
  public:
    BS(double sigma, double r, const double K, const double T, const int N_t, const int N_S);

    void step(int N);     // n time steps of backward evolution
    void step_linearCC(int N); // evolution with linear boundary conditions
    void step_quadraticCC(int N); // evolution with quadratic boundary conditions

    vector<double> get_C() const {return C;}
    vector<double> get_C_analitic(); // returns the analytic value of the call option at time t
    double get_error_L1(); // see Anwar's article for definition, uses L1 norm
    double get_error_inf(); // error with infinity norm (takes only the maximum error in S)
    vector<double> get_error(); // vector of deviations from the analytic solution
    double get_t() const {return t;}
    void save_data(string file_name) const;     // methods to save data in the Data directory
    void save_grid(string file_name) const;
    void save_analitic(string file_name);
    void save_error(vector<double> error, string file_name);  // saves the error vector (to be calculated)
    void save(string file_name) const {save_data(file_name); save_grid(file_name);}

  private:
    double sigma;   // volatility (may depend on t)
    double r;       // riskless rate (may depend on t)
    const double K;       // strike price
    const double T;       // maturity
    const int N_t;        // temporal discretization steps
    const int N_S;        // discretization steps for the asset value
    const double S_max;   // maximum asset value, set to boundary * K
    vector<double> S;       // grid of returns
    double t;             // current time, initialized to T (backward evolution)
    const double delta_t;
    const double delta_S;
    vector<double> C;    // call option value at time t for each S
    double C_aux;        // call option value at S_tilda / 2 (for interpolation)
    double linear_boundary(); // approximates the function linearly beyond the boundary
    double quadratic_boundary();  // approximates the function with a quadratic polynomial beyond the boundary,
                                  // in both cases interpolation is done with a point at S_tilda
    void smoothing();    // function that smooths the corner point with an exponential

    double c1(int i);     // coefficient of C(i, t)
    double c2(int i);     // coefficient of C(i+1, t)
    double c3(int i);     // coefficient of C(i-1, t)
};

double power(double b, int e);    // = b^e

void print(vector<double> v);

double normal_cum(double x);  // cumulative normal function with error O(10^-15) = double precision error

double normal(double x); // normal function

double simpson_int(vector<double> f, const double dx); // f vector containing f evaluated on the grid

bool iseven(int n); // true if n is even
