#include "ResultStepOne.h"
#include "seasonality.h"
#include <math.h>
#include <vector>
using namespace std;

double*** phi_hat_element(int Nstep, int Nitera, int h2, double lambda, double dt, int A, int C, int K, double pi, double p1, double f1, double** J_tree) {
	double*** phi_hat_tree = new double** [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		phi_hat_tree[i] = new double* [Nstep + 1];
		for (int j = 0; j < Nstep + 1; j++) {
			phi_hat_tree[i][j] = new double[Nitera];
		}
	}
	for (int i = 0; i < Nstep + 1; i++)
		for (int j = 0; j < Nitera; j++) {
			phi_hat_tree[Nstep][i][j] = h2;
		}
	for (int i = 0; i < Nstep; i++)
		for (int j = 0; j < i + 1; j++) {
			phi_hat_tree[i][j][Nitera - 1] = 0;
		}
	double** cond_phi_hat_tree = new double* [Nstep];
	for (int i = 0; i < Nstep; i++) {
		cond_phi_hat_tree[i] = new double[Nstep];
	}

	for (int i = Nstep - 1; i > -1; i--)
		for (int j = 0; j < i + 1; j++)
			for (int k = 1; k < Nitera; k++) {
				cond_phi_hat_tree[i][j] = (1 - exp(-lambda / Nstep)) * phi_hat_tree[i + 1][i + 1][k] + exp(-lambda / Nstep) * phi_hat_tree[i + 1][j][k];
				phi_hat_tree[i][j][k] = (cond_phi_hat_tree[i][j] + dt * (C + pow(phi_hat_tree[i][j][k - 1], 2) / (A + K + (1 - pi) * p1 + f1 * J_tree[i][j]))) / (1 + (dt * 2 * phi_hat_tree[i][j][k - 1] / (A + K + (1 - pi) * p1 + f1 * J_tree[i][j])));
			}
	return phi_hat_tree;
}

double*** psi_hat_element(int Nstep, int h1, double lambda, double A, double p0, double p1, int f0, int f1, int K, double pi, int Nitera, double dt, double bar_Alpha, double** bar_Q_tree, double** Q_r_tree, double** J_tree, double*** phi_hat_tree, vector<double> season) {
	double*** psi_hat_tree = new double** [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		psi_hat_tree[i] = new double* [Nstep + 1];
		for (int j = 0; j < Nstep + 1; j++) {
			psi_hat_tree[i][j] = new double[Nstep + 1];
		}
	}
	for (int i = 0; i < Nstep + 1; i++)
		for (int j = 0; j < Nstep + 1; j++) {
			psi_hat_tree[Nstep][i][j] = h1;
		}
	double*** cond_psi_hat_tree = new double** [Nstep];
	for (int i = 0; i < Nstep + 1; i++) {
		cond_psi_hat_tree[i] = new double* [Nstep];
		for (int j = 0; j < Nstep + 1; j++) {
			cond_psi_hat_tree[i][j] = new double[Nstep];
		}
	}
	int bar_P;
	for (int i = Nstep - 1; i > -1; i--)
		for (int j = 0; j < i + 1; j++)
			for (int k = 0; k < i + 1; k++) {
				cond_psi_hat_tree[i][j][k] = 0.5 * (1 - exp(-lambda / Nstep)) * psi_hat_tree[i + 1][j][i + 1] + 0.5 * exp(-lambda / Nstep) * psi_hat_tree[i + 1][j][k] + 0.5 * (1 - exp(-lambda / Nstep)) * psi_hat_tree[i + 1][j + 1][i + 1] + 0.5 * exp(-lambda / Nstep) * psi_hat_tree[i + 1][j + 1][k];
				bar_P = p0 + K * bar_Q_tree[i][j] + p1 * ((1 - pi) * bar_Q_tree[i][j] + pi * Q_r_tree[i][j]) + (f0 + f1 * (bar_Q_tree[i][j] - season[i] - bar_Alpha)) * J_tree[i][k];
				psi_hat_tree[i][j][k] = (cond_psi_hat_tree[i][j][k] - dt * phi_hat_tree[i][k][Nitera - 1] * bar_P / (A + K + (1 - pi) * p1 + f1 * J_tree[i][k])) / (1 + (dt * phi_hat_tree[i][k][Nitera - 1] / (A + K + (1 - pi) * p1 + f1 * J_tree[i][k])));
			}
	return psi_hat_tree;
}


vector<double> phi_hat(int n, double*** phi_hat_tree, double lambda, int h2, int Nitera, vector<double> randnPoi) {
	vector<double> phi_hat_simulation(n + 1, 0.);
	phi_hat_simulation[n] = h2;
	phi_hat_simulation[0] = phi_hat_tree[0][0][Nitera - 1];
	int indexJ = 0;
	for (int i = 1; i < n; i++) {
		if (randnPoi[i - 1] == exp(-lambda / n)) indexJ = i;
		phi_hat_simulation[i] = phi_hat_tree[i][indexJ][Nitera - 1];
	}
	return phi_hat_simulation;
}

vector<double> psi_hat(int n, double*** psi_hat_tree, double lambda, int h1, int Nitera, vector<double> randnPoi, vector<double> randnCom) {
	vector<double> psi_hat_simulation(n + 1, 0.);
	psi_hat_simulation[n] = h1;
	psi_hat_simulation[0] = psi_hat_tree[0][0][0];
	int indexCom = 0;
	int indexJ = 0;
	for (int i = 1; i < n; i++) {
		if (randnPoi[i - 1] == exp(-lambda / n)) indexJ = i;
		if (randnCom[i - 1] == 1) indexCom += 1;
		psi_hat_simulation[i] = psi_hat_tree[i][indexCom][indexJ];
	}
	return psi_hat_simulation;
}

vector<double> bar_P(int Nstep, double p0, double p1, int K, double pi, vector<double> bar_Q_simulation, vector<double> Q_r_simulation, double f0, double f1, double bar_Alpha, double dt, vector<double> J_simulation, vector<double> season) {
	vector<double> bar_P_simulation(Nstep + 1, 0.);
	for (int i = 0; i < Nstep + 1; i++) {
		bar_P_simulation[i] = p0 + K * bar_Q_simulation[i] + p1 * ((1 - pi) * bar_Q_simulation[i] + pi * Q_r_simulation[i]) + (f0 + f1 * (bar_Q_simulation[i] - season[i] - bar_Alpha)) * J_simulation[i];
	}
	return bar_P_simulation;
}


vector<double> S_hat(int n, vector<double> phi_hat_simulation, vector<double> psi_hat_simulation, vector<double> bar_P_simulation, vector<double> J_simulation, int A, double p1, int f1, double dt, double S_init, double pi, int K) {
	vector<double> TildeA_simulation(n + 1, 0.);
	for (int i = 0; i < n + 1; i++) {
		TildeA_simulation[i] = -(psi_hat_simulation[i] + bar_P_simulation[i]) / (A + K + (1 - pi) * p1 + f1 * J_simulation[i]);
	}
	vector<double> integral1(n + 1, 0.);
	for (int i = 1; i < n + 1; i++) {
		integral1[i] = integral1[i - 1] + phi_hat_simulation[i] / (A + K + (1 - pi) * p1 + f1 * J_simulation[i]) * dt;
	}
	vector<double> integral2(n + 1, 0.);
	for (int i = 1; i < n + 1; i++) {
		integral2[i] = integral2[i - 1] + TildeA_simulation[i] * exp(integral1[i]) * dt;
	}
	vector<double> S_hat_simulation(n + 1, 0.);
	S_hat_simulation[0] = S_init;
	for (int i = 1; i < n + 1; i++) {
		S_hat_simulation[i] = exp(-integral1[i]) * (S_init + integral2[i]);
	}
	return S_hat_simulation;
}





vector<double> alpha_hat(int n, vector<double> phi_hat_simulation, vector<double> S_hat_simulation, vector<double> psi_hat_simulation, vector<double> bar_P_simulation, vector<double> J_simulation, int A, int p1, int f1, double pi, int K) {
	vector<double> alpha_hat_simulation(n + 1, 0.);
	for (int i = 0; i < n + 1; i++) {
		alpha_hat_simulation[i] = -(psi_hat_simulation[i] + phi_hat_simulation[i] * S_hat_simulation[i] + bar_P_simulation[i]) / (A + K + (1 - pi) * p1 + f1 * J_simulation[i]);
	}
	return alpha_hat_simulation;
}


