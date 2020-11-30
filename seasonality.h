#pragma once

#include <vector>
using namespace std;


/*the trajectory of Q_r*/
vector<double> Q_r(int n, vector<double> path, double dt, double sigma, vector<double> season);

/*the trajectory of bar_Q*/
vector<double> bar_Q(int n, vector<double> path, double dt, double sigma, vector<double> season);

/*the trajectory of Q*/
vector<double> Q(int n, vector<double> pathCom, vector<double> pathInd, double dt, double sigmaCom, double sigmaInd, vector<double> season);

/*the trajectory of Y*/
vector<double> Y(int n, double delay, vector<double> pathPoi, double dt, double lambda);

/*the trajectory of J*/
vector<double> J(int n, vector<double> pathY, double delay);

/*the tree of Q_r*/
double** Q_r_element(int Nstep, double sigma, double dt, vector<double> season);

/*the tree of bar_Q*/
double** bar_Q_element(int Nstep, double sigma, double dt, vector<double> season);

/*the tree of Q*/
double*** Q_element(int Nstep, double sigmaCom, double sigmaInd, double dt, vector<double> season);

/*the tree of J*/
double** J_element(int Nstep, double dt, double delay);