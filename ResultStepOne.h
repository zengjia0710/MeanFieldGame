#pragma once

#include <vector>
using namespace std;

/*construct the tree of phi_hat*/
double*** phi_hat_element(int Nstep, int Nitera, int h2, double lambda, double dt, int A, int C, int K, double pi, double p1, double f1, double** J_tree);

/*construct the tree of psi_hat*/
double*** psi_hat_element(int Nstep, int h1, double lambda, double A, double p0, double p1, int f0, int f1, int K, double pi, int Nitera, double dt, double bar_Alpha, double** bar_Q_tree, double** Q_r_tree, double** J_tree, double*** phi_hat_tree, vector<double> season);

/*pick the trajectory of phi_hat*/
vector<double> phi_hat(int n, double*** phi_hat_tree, double lambda, int h2, int Nitera, vector<double> randnPoi);

/*pick the trajectory of psi_hat*/
vector<double> psi_hat(int n, double*** psi_hat_tree, double lambda, int h1, int Nitera, vector<double> randnPoi, vector<double> randnCom);

/*pick the trajectory of bar_P*/
vector<double> bar_P(int Nstep, double p0, double p1, int K, double pi, vector<double> bar_Q_simulation, vector<double> Q_r_simulation, double f0, double f1, double bar_Alpha, double dt, vector<double> J_simulation, vector<double> season);

/*pick the trajectory of S_hat*/
vector<double> S_hat(int n, vector<double> phi_hat_simulation, vector<double> psi_hat_simulation, vector<double> bar_P_simulation, vector<double> J_simulation, int A, double p1, int f1, double dt, double S_init, double pi, int K);

/*pick the trajectory of alpha_hat*/
vector<double> alpha_hat(int n, vector<double> phi_hat_simulation, vector<double> S_hat_simulation, vector<double> psi_hat_simulation, vector<double> bar_P_simulation, vector<double> J_simulation, int A, int p1, int f1, double pi, int K);