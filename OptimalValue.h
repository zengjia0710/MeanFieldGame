#pragma once

#include <vector>
#include "seasonality.h"
#include "rand.h"
using namespace std;

double OptimalValue(int Nstep, int Nmc, int Nitera, vector<double> season, vector<double> phi_simulation, int A, int C, int K, double p0, double p1, double bar_Alpha, int f0, int f1, int h0, int h1, int h2, double dt, double pi, double S_init, double delay, double lambda, double sigmaCom, double sigmaInd, double*** phi_hat_tree, double*** psi_hat_tree);
