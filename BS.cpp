
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>
#include "BS.h"
using namespace std;

double normal_cum(double x) {
  const double x0 = -8;
  if(x <= x0) return 0;
  if(x > 8) return 1;
  const double dx = pow(10, -4);
  long int N = round( (x - x0)/dx );
  if( !iseven(N) ) N = N - 1;  // N must be even so that N+1 is odd (see line below)
  vector<double> f(N+1);

  for(int i = 0; i < f.size(); i++) f[i] = normal( x0 + i*dx );

  return simpson_int(f, dx);
}

double normal(double x) {
  return exp(- x*x / 2)/sqrt(2 * M_PI);
}

double simpson_int(vector<double> f, const double dx) {

  for(int i = 1; i < f.size() - 1; i++) {
    if(iseven(i)) f[i] *= 2;
    else f[i] *= 4;
  }

  double sum = 0;
  sum = accumulate(f.begin(), f.end(), sum);
  sum *= dx/3;
  return sum;
}

bool iseven(int n) {
  if( n%2 == 0) return true;
  else return false;
}

BS::BS(double sigmasigma, double rr, double KK, double TT, int N_tt, int N_SS)
: sigma(sigmasigma), r(rr), K(KK), T(TT), N_t(N_tt), N_S(N_SS), S_max(boundary*K), t(T), delta_t(T/N_t), delta_S(S_max/N_S)
{
  C.resize(N_S);
  for(int i = 0; i < N_S; i++)      // boundary conditions
  {
    C[i] = max(i*delta_S-K, 0.0);
  }

  C_aux = S_tilda/2 - K;          // useful only for step_quadraticCC

  S.resize(N_S);
  for(int i = 0; i < N_S; i++)
  {
    S[i] = i * delta_S;
  }
}

void BS::smoothing() {
  double delta = 5;

  for(int i = 0; i < C.size(); i++) {
    if( S[i] >= (K - delta) && S[i] <= K ) {
      C[i] = exp( ( S[i] - K - 1 ) );
    }
  }
}


void BS::step(int N)
{
  vector<double> aux(N_S);    // to avoid overwriting C
  aux[N_S-1] = C[N_S-1];      // boundary condition, aux[0] is already = 0
  for(int n = 0; n < N; n++) {
    t -= delta_t;
    for(int i = 1; i < N_S - 1; i++) aux[i] = c1(i) * C[i+1] + c2(i) * C[i] + c3(i) * C[i-1];
    C = aux;
  }
}


void BS::step_linearCC(int N)
{
  vector<double> aux(N_S);    // to avoid overwriting C
  for(int n = 0; n < N; n++) {
    t -= delta_t;
    for(int i = 1; i < N_S - 1; i++) aux[i] = c1(i) * C[i+1] + c2(i) * C[i] + c3(i) * C[i-1];
    aux[N_S-1] = linear_boundary();
    C = aux;
  }
}

void BS::step_quadraticCC(int N)
{
  vector<double> aux(N_S);    // to avoid overwriting C
  for(int n = 0; n < N; n++) {
    t -= delta_t;
    for(int i = 1; i < N_S - 1; i++) aux[i] = c1(i) * C[i+1] + c2(i) * C[i] + c3(i) * C[i-1];
    aux[N_S-1] = quadratic_boundary();
    C = aux;
  }
}


double BS::quadratic_boundary() {
  double c1; double c2; double c3;
  double phi; double psi; double xi;
  double Aq; double Bq; double Cq;

  c1 = (S[N_S-2] - S_tilda/2) * (S[N_S-2] - S_tilda);
  c2 = -S_tilda / 2 * (S_tilda/2 - S[N_S-2]);
  c3 = S_tilda / 2 * (S_tilda - S[N_S-2]);

  phi = C[N_S-2] / c1;
  psi = C_aux / c2;
  xi = (S_tilda - K) / c3;

  Aq = phi + psi + xi;
  Bq = - ( 3/2 * phi * S_tilda + psi * (S_tilda + S[N_S-2]) + xi * (S[N_S-2] + S_tilda / 2) );
  Cq = phi / 2 * pow(S_tilda, 2) + psi * S_tilda * S[N_S-2] + xi * S_tilda / 2 * S[N_S-2];

  C_aux = (1 - delta_t*r) * C_aux
          + delta_t/2 * ( Aq*(pow(sigma, 2)/2 + r)*pow(S_tilda, 2) + Bq*r*S_tilda);

  return Aq * pow(S[N_S-1], 2) + Bq * S[N_S-1] + Cq;

}

double BS::linear_boundary() {
  double bound = (S_tilda - K - C[N_S-2]) / (S_tilda - S[N_S-2]) * delta_S + C[N_S-2];
  return bound;
}

void BS::save_data(string file_name) const {
  fstream file;
  file.open("Data/" + file_name, ios_base::out);
  for(int i = 0; i < N_S * see_till / boundary; i++) file << C[i] << "\t";
  file.close();
}

void BS::save_grid(string file_name) const {
  fstream file;
  file.open("Data/grid_" + file_name, ios_base::out);
  for(int i = 0; i < N_S * see_till / boundary; i++) file << S[i] << "\t";
  file.close();
}

void BS::save_analitic(string file_name) {
  vector<double> C_analitic(N_S);
  C_analitic = get_C_analitic();
  fstream file;
  file.open("Data/analitic_" + file_name, ios_base::out);
  for(int i = 0; i < N_S * see_till / boundary; i++) file << C_analitic[i] << "\t";
  file.close();
}

void BS::save_error(vector<double> error, string file_name) {
  fstream file;
  file.open("Data/" + file_name, ios_base::out);
  for(int i = 0; i < error.size(); i++) file << error[i] << "\t";
  file.close();
}

vector<double> BS::get_C_analitic() {
  double d1;
  double d2;
  vector<double> C_analitic(N_S * see_till / boundary);
  for(int i = 0; i < C_analitic.size(); i++) {
    if( S[i] == 0 ) d1 = -8; // avoid log divergences, x = -8 = -inf for normal_cum
    else d1 = ( log(S[i]/K) + (r  + pow(sigma, 2) / 2) * (T - t) ) / ( sigma * sqrt(T - t) );
    d2 = d1 - sigma * sqrt(T - t);
    C_analitic[i] = S[i] * normal_cum(d1) - K * exp(- r * (T - t) ) * normal_cum(d2);
  }
  return C_analitic;
}

double BS::get_error_L1() {
  vector<double> C_analitic;
  C_analitic = get_C_analitic();

  double C_an_norm1 = 0;
  C_an_norm1 = accumulate(C_analitic.begin(), C_analitic.end(), C_an_norm1);

  vector<double> error = get_error();

  double error_L1 = 0;
  error_L1 = accumulate(error.begin(), error.end(), error_L1);
  error_L1 /= C_an_norm1;

  return error_L1;
}

double BS::get_error_inf() {
  vector<double> error = get_error();
  double error_inf = *max_element(error.begin(), error.end());

  return error_inf;
}

vector<double> BS::get_error() {
  vector<double> error(N_S * see_till / boundary);
  error = get_C_analitic();
  for(int i = 0; i < error.size(); i++) {
    error[i] = abs(error[i] - C[i]);
  }

  return error;
}

double BS::c1(int i)
{
  return (r * i + sigma*sigma * i*i ) * delta_t / 2;
}

double BS::c2(int i)
{
  return 1 - delta_t * r - sigma*sigma * delta_t * i*i;
}

double BS::c3(int i)
{
  return (sigma*sigma * i*i - r * i) * delta_t / 2;
}

double power(double b, int e)
{
  double pow = 1;
  for(int i = 0; i<e; i++)
  {
    pow *= b;
  }
  return pow;
}

void print(vector<double> v)
{
  cout << "(";
  for(int i = 0; i<v.size(); i++) {
    if(i!=v.size()-1) cout << v[i] << ", ";
    else cout << v[i] << ")" << endl;
  }
}
