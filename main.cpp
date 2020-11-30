#include <iostream>
#include <Python.h>
#include "plot.h"
#include "seasonality.h"
#include "rand.h"
#include "ResultStepOne.h"
#include "ResultStepTwo.h"
#include "OptimalValue.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <ctime>
using namespace std;

/*parameters*/
double T = 2;
const int Nstep = 95;
double dt = T / Nstep;
int Nitera = 100;
int Nmc = 5000;

int A = 150;
int C = 80;
int f0 = 0;
int f1 = 10000;
double p0 = 6.159423723;
double p1 = 87.4286117;
int h0 = 0;
int h1 = 0;
int h2 = 100;
double bar_Alpha = -0.1;
double delay = 0.12;
double sigma[] = { 0.31, 0.56 };
double lambda = 2;
double pi = 0.5;
int K = 50;
double S_init = -0.5;

int main() {
	/*seasonality £¨from real data)*/
	double data[] = { 0.26759617, 0.24771933, 0.23588383, 0.221369, 0.21174, 0.2047625, 0.20651067, 0.20098083, 0.20826067, 0.22095067,
	   0.24346833, 0.27283267, 0.3382265, 0.42920433, 0.4875495, 0.50948433, 0.487712, 0.4537295, 0.40911717, 0.3728925,
	   0.347346, 0.3419715, 0.32684, 0.320009, 0.32065767, 0.32586567, 0.31492483, 0.31607417, 0.30411783, 0.29950567,
	   0.307519, 0.33259367, 0.375465, 0.45608333, 0.599178,0.70970583, 0.7364855, 0.736731, 0.70612667, 0.67284583,
	   0.66692767, 0.64925583, 0.604485, 0.55684567, 0.515597, 0.45097333, 0.3822625, 0.31841833 };

	vector<double> season(Nstep + 1, 0.);
	for (int i = 0; i < 48; i++) {
		season[i] = data[i];
	}
	// if T > 1
	if (Nstep / 47 > 1) {
		for (int i = 48; i < Nstep + 1; i++) {
			season[i] = data[i % 48];
		}
	}


	/*construct the tree of Q_r*/
	double** Q_r_tree = new double* [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		Q_r_tree[i] = new double[Nstep + 1];
	}
	Q_r_tree = Q_r_element(Nstep, sigma[0], dt, season);


	/*construct the tree of bar_Q*/
	double** bar_Q_tree = new double* [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		bar_Q_tree[i] = new double[Nstep + 1];
	}
	bar_Q_tree = bar_Q_element(Nstep, sigma[0], dt, season);

	/*construct the tree of Q*/
	double*** Q_tree = new double** [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		Q_tree[i] = new double* [Nstep + 1];
		for (int j = 0; j < Nstep + 1; j++) {
			Q_tree[i][j] = new double[Nstep + 1];
		}
	}
	Q_tree = Q_element(Nstep, sigma[0], sigma[1], dt, season);

	/*construct the tree of J*/
	double** J_tree = new double* [Nstep + 1];
	J_tree = J_element(Nstep, dt, delay);


	/*construct the tree of phi_hat*/
	double*** phi_hat_tree = new double** [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		phi_hat_tree[i] = new double* [Nstep + 1];
		for (int j = 0; j < Nstep + 1; j++) {
			phi_hat_tree[i][j] = new double[Nitera];
		}
	}
	phi_hat_tree = phi_hat_element(Nstep, Nitera, h2, lambda, dt, A, C, K, pi, p1, f1, J_tree);

	/*construct the tree of psi_hat*/
	double*** psi_hat_tree = new double** [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		psi_hat_tree[i] = new double* [Nstep + 1];
		for (int j = 0; j < Nstep + 1; j++) {
			psi_hat_tree[i][j] = new double[Nstep + 1];
		}
	}
	psi_hat_tree = psi_hat_element(Nstep, h1, lambda, A, p0, p1, f0, f1, K, pi, Nitera, dt, bar_Alpha, bar_Q_tree, Q_r_tree, J_tree, phi_hat_tree, season);


	/*generate the random path of common noise*/
	vector<double> randnCom(Nstep, 0.);
	randnCom = randBM(Nstep);


	/*generate the random path of independent noise*/
	vector<double> randnInd(Nstep, 0.);
	randnInd = randBM(Nstep);

	/*generate the random path of Poisson process*/
	vector<double> randnPoi(Nstep, 0.);
	randnPoi = randPoisson(Nstep, lambda, Nstep);


	/*the trajectories of Q_r, bar_Q, Q, Y and J*/
	vector<double> Q_r_simulation(Nstep + 1, 0.);
	Q_r_simulation = Q_r(Nstep, randnCom, dt, sigma[0], season);

	vector<double> bar_Q_simulation(Nstep + 1, 0.);
	bar_Q_simulation = bar_Q(Nstep, randnCom, dt, sigma[0], season);

	vector<double> Q_simulation(Nstep + 1, 0.); //Agent 1
	Q_simulation = Q(Nstep, randnCom, randnInd, dt, sigma[0], sigma[1], season);

	vector<double> Y_simulation(Nstep + 1, 0.);
	Y_simulation = Y(Nstep, delay, randnPoi, dt, lambda);

	vector<double> J_simulation(Nstep + 1, 0.);
	J_simulation = J(Nstep, Y_simulation, delay);

	clock_t time_start = clock();
	/*the trajectory of phi_hat*/
	vector<double> phi_hat_simulation(Nstep + 1, 0.);
	phi_hat_simulation = phi_hat(Nstep, phi_hat_tree, lambda, h2, Nitera, randnPoi);

	/*the trajectory of psi_hat*/
	vector<double> psi_hat_simulation(Nstep + 1, 0.);
	psi_hat_simulation = psi_hat(Nstep, psi_hat_tree, lambda, h1, Nitera, randnPoi, randnCom);

	/*the trajectory of bar_P*/
	vector<double> bar_P_simulation(Nstep + 1, 0.);
	bar_P_simulation = bar_P(Nstep, p0, p1, K, pi, bar_Q_simulation, Q_r_simulation, f0, f1, bar_Alpha, dt, J_simulation, season);

	/*the trajectory of S_hat*/
	vector<double> S_hat_simulation(Nstep + 1, 0.);
	S_hat_simulation = S_hat(Nstep, phi_hat_simulation, psi_hat_simulation, bar_P_simulation, J_simulation, A, p1, f1, dt, S_init, pi, K);

	/*the trajectory of alpha_hat*/
	vector<double> alpha_hat_simulation(Nstep + 1, 0.);
	alpha_hat_simulation = alpha_hat(Nstep, phi_hat_simulation, S_hat_simulation, psi_hat_simulation, bar_P_simulation, J_simulation, A, p1, f1, pi, K);

	clock_t time_end = clock();
	cout << "time use of psi_hat:" << 1000 * (double)(time_end - time_start) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << "phi_hat[0] = " << phi_hat_simulation[0] << endl;
	cout << "psi_hat[0] = " << psi_hat_simulation[0] << endl;
	cout << "alpha_hat[0] = " << alpha_hat_simulation[0] << endl;


	/*Price*/
	vector<double> price(Nstep + 1, 0.);
	for (int i = 0; i < Nstep + 1; i++) {
		price[i] = p0 + 87.4286117 * pi * Q_r_simulation[i] + 87.4286117 * (1 - pi) * (bar_Q_simulation[i] + alpha_hat_simulation[i]);
	}


	/*the trajectory of phi*/
	vector<double> phi_simulation(Nstep + 1, 0.);
	phi_simulation = phi(Nstep, h2, C, A, dt);

	/*the trajectory of phi*/
	clock_t time_start_psi = clock();
	vector<double> psi_simulation(Nstep + 1, 0.); //Agent 1
	psi_simulation[Nstep] = h1;
	for (int i = 0; i < Nstep; i++) {
		psi_simulation[i] = psi(i, Nstep, Nitera, Nmc, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, randnCom, randnInd, randnPoi, dt, sigma[0], sigma[1], phi_hat_tree, psi_hat_tree, phi_simulation, season);
	}


	/*the trajectory of P*/
	vector<double> P_simulation(Nstep + 1, 0.); //Agent 1
	P_simulation = P_StepTwo(Nstep, Nitera, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, randnCom, randnInd, randnPoi, dt, sigma[0], sigma[1], phi_hat_tree, psi_hat_tree, season);

	/*the trajectory of S*/
	vector<double> S_simulation(Nstep + 1, 0.); //Agent 1
	S_simulation = S(Nstep, phi_simulation, psi_simulation, P_simulation, A, K, dt, S_init);

	/*the trajectory of alpha*/
	vector<double> alpha_simulation(Nstep + 1, 0.); //Agent 1
	alpha_simulation = alpha(Nstep, phi_simulation, psi_simulation, S_simulation, P_simulation, A, K);

	clock_t time_end_psi = clock();
	cout << "time use of psi:" << 1000 * (double)(time_end_psi - time_start_psi) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << "psi[0] = " << psi_simulation[0] << endl;
	cout << "alpha[0] = " << alpha_simulation[0] << endl;
	system("pause");


	/*Plotting*/

	double a[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { a[i] = i; }
	double b[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { b[i] = bar_Q_simulation[i]; }
	double c[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { c[i] = Q_simulation[i]; }
	double d[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { d[i] = J_simulation[i]; }
	double e[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { e[i] = bar_Q_simulation[i] + alpha_hat_simulation[i]; }
	double e1[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { e1[i] = Q_simulation[i] + alpha_simulation[i]; }
	double f1[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { f1[i] = psi_simulation[i]; }
	double g[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { g[i] = alpha_hat_simulation[i]; }
	double g1[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { g1[i] = alpha_simulation[i]; }
	double h[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { h[i] = S_hat_simulation[i]; }
	double h1[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { h1[i] = S_simulation[i]; }
	double price1[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { price1[i] = p0 + 87.4286117 * bar_Q_simulation[i]; }
	double price2[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { price2[i] = p0 + 87.4286117 * (pi * Q_r_simulation[i] + (1 - pi) * (bar_Q_simulation[i] + alpha_hat_simulation[i])); }
	double Y_hat[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { Y_hat[i] = phi_hat_simulation[i] * S_hat_simulation[i] + psi_hat_simulation[i]; }
	double Y[Nstep + 1] = {};
	for (int i = 0; i < Nstep + 1; ++i) { Y[i] = phi_simulation[i] * S_simulation[i] + psi_simulation[i]; }
	pythonInitial();
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.argv = ['python.py']");
	plot(a, Nstep + 1, b, "hat_Q", "trajectories of hat_Q and Q", c, "Q");
	plot(a, Nstep + 1, "trajectories of J", d);
	plot(a, Nstep + 1, e, "hat_Q + hat_alpha", "trajectories of Q + alpha and hat_Q + hat_alpha", e1, "Q + alpha");
	plot(a, Nstep + 1, h, "hat_S", "trajectories of hat_S and S", h1, "S");
	plot(a, Nstep + 1, g, "hat_alpha", "trajectories of hat_alpha and alpha", g1, "alpha");
	plot(a, Nstep + 1, price1, "price without alpha", "prices", price2, "price with alpha (pi = 0.5)");
	Py_Finalize();
	system("pause");
}