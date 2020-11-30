#pragma once

#include <vector>
using namespace std;
#include "rand.h"

/*pick the trajectory of P*/
vector<double> P_StepTwo(int Nstep, int Nitera, int A, int K, double p0, double p1, double pi, int f0, int f1, double h1, double h2, double S_init, double bar_Alpha, double delay, double lambda, vector<double> pathCom, vector<double> pathInd, vector<double> pathJump, double dt, double sigmaCom, double sigmaInd, double*** phi_hat_tree, double*** psi_hat_tree, vector<double> season);

/*pick the trajectory of phi*/
vector<double> phi(int Nstep, int h2, int C, int A, double dt);

/*compute the value of psi at time k*/
double psi(int timeK, int Nstep, int Nitera, int Nmc, int A, int K, double p0, double p1, double pi, int f0, int f1, double h1, double h2, double S_init, double bar_Alpha, double delay, double lambda, vector<double> pathCom, vector<double> pathInd, vector<double> pathJump, double dt, double sigmaCom, double sigmaInd, double*** phi_hat_tree, double*** psi_hat_tree, vector<double> phi_simulation, vector<double> season);

/*pick the trajectory of S*/
vector<double> S(int Nstep, vector<double> phi_simulation, vector<double> psi_simulation, vector<double> P_simulation, int A, int K, double dt, double S_init);

/*pick the trajectory of alpha*/
vector<double> alpha(int Nstep, vector<double> phi_simulation, vector<double> psi_simulation, vector<double> S_simulation, vector<double> P_simulation, int A, int K);
