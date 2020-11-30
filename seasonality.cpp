#include "seasonality.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;


vector<double> Q_r(int Nstep, vector<double> path, double dt, double sigma, vector<double> season) {
	vector<double> ret(Nstep + 1, 0.);
	ret[0] = season[0];
	int sum_randn = 0;
	double t = 0;
	for (int i = 1; i < Nstep + 1; i++) {
		sum_randn += path[i - 1];
		t += dt;
		ret[i] = season[i] * exp(-pow(sigma, 2) / 2 * t) * exp(sigma * sqrt(dt) * sum_randn);
	}
	return ret;
}

vector<double> bar_Q(int Nstep, vector<double> path, double dt, double sigma, vector<double> season) {
	vector<double> ret(Nstep + 1, 0.);
	ret[0] = season[0];
	int sum_randn = 0;
	double t = 0;
	for (int i = 1; i < Nstep + 1; i++) {
		sum_randn += path[i - 1];
		t += dt;
		ret[i] = season[i] * exp(-pow(sigma, 2) / 2 * t) * exp(sigma * sqrt(dt) * sum_randn);
	}
	return ret;
}

vector<double> Q(int Nstep, vector<double> pathCom, vector<double> pathInd, double dt, double sigmaCom, double sigmaInd, vector<double> season) {
	vector<double> ret(Nstep + 1, 0.);
	ret[0] = season[0];
	int sum_randn1 = 0;
	int sum_randn2 = 0;
	double t = 0;
	for (int i = 1; i < Nstep + 1; i++) {
		sum_randn1 += pathCom[i - 1];
		sum_randn2 += pathInd[i - 1];
		t += dt;
		ret[i] = season[i] * exp((-pow(sigmaCom, 2) / 2 - pow(sigmaInd, 2) / 2) * t) * exp(sigmaCom * sqrt(dt) * sum_randn1) * exp(sigmaInd * sqrt(dt) * sum_randn2);
	}
	return ret;
}


vector<double> Y(int Nstep, double delay, vector<double> pathPoi, double dt, double lambda) {
	vector<double> ret(Nstep + 1, 0.);
	ret[0] = 2 * delay;
	for (int i = 1; i < Nstep + 1; i++) {
		ret[i] = ret[i - 1] + dt * (1 - lambda * ret[i - 1]) - ret[i - 1] * pathPoi[i - 1];
	}
	return ret;
}

vector<double> J(int Nstep, vector<double> pathY, double delay) {
	vector<double> ret(Nstep + 1, 0.);
	for (int i = 0; i < Nstep + 1; i++) {
		if (pathY[i] < delay) ret[i] = 1;
		else ret[i] = exp(-pow(10, 9) * pow(pathY[i] - delay, 2));
	}
	return ret;
}

double** Q_r_element(int Nstep, double sigma, double dt, vector<double> season) {
	double** Q_r_tree = new double* [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		Q_r_tree[i] = new double[Nstep + 1];
	}
	Q_r_tree[0][0] = season[0];
	for (int i = 0; i < Nstep; i++)
		for (int j = 0; j < i + 2; j++) {
			double t = (i + 1) * dt;
			Q_r_tree[i + 1][j] = season[i + 1] * exp(-pow(sigma, 2) / 2 * t) * exp(sigma * sqrt(dt) * (-i - 1 + 2 * j));
		}
	return Q_r_tree;
}

double** bar_Q_element(int Nstep, double sigma, double dt, vector<double> season) {
	double** bar_Q_tree = new double* [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		bar_Q_tree[i] = new double[Nstep + 1];
	}
	bar_Q_tree[0][0] = season[0];
	for (int i = 0; i < Nstep; i++)
		for (int j = 0; j < i + 2; j++) {
			double t = (i + 1) * dt;
			bar_Q_tree[i + 1][j] = season[i + 1] * exp(-pow(sigma, 2) / 2 * t) * exp(sigma * sqrt(dt) * (-i - 1 + 2 * j));
		}
	return bar_Q_tree;
}

double*** Q_element(int Nstep, double sigmaCom, double sigmaInd, double dt, vector<double> season) {
	double*** Q_tree = new double** [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		Q_tree[i] = new double* [Nstep + 1];
		for (int j = 0; j < Nstep + 1; j++) {
			Q_tree[i][j] = new double[Nstep + 1];
		}
	}
	Q_tree[0][0][0] = season[0];
	for (int i = 0; i < Nstep; i++)
		for (int j = 0; j < i + 2; j++)
			for (int k = 0; k < i + 2; k++) {
				double t = (i + 1) * dt;
				Q_tree[i + 1][j][k] = season[i + 1] * exp((-pow(sigmaCom, 2) / 2 - pow(sigmaInd, 2) / 2) * t) * exp(sigmaCom * sqrt(dt) * (-i - 1 + 2 * j)) * exp(sigmaInd * sqrt(dt) * (-i - 1 + 2 * k));
			}
	return Q_tree;
}

double** J_element(int Nstep, double dt, double delay) {
	/*construct the tree of Y*/
	double** Y_tree = new double* [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		Y_tree[i] = new double[Nstep + 1];
	}
	for (int i = 0; i < Nstep + 1; i++) {
		Y_tree[i][0] = 2 * delay + i * dt;
	}
	for (int i = 2; i < Nstep + 1; i++)
		for (int j = 1; j < i; j++) {
			Y_tree[i][j] = (i - j) * dt;
		}

	/*construct the tree of J*/
	double** J_tree = new double* [Nstep + 1];
	for (int i = 0; i < Nstep + 1; i++) {
		J_tree[i] = new double[Nstep + 1];
	}
	for (int i = 0; i < Nstep + 1; i++)
		for (int j = 0; j < i + 1; j++) {
			if (Y_tree[i][j] < delay) J_tree[i][j] = 1;
			else J_tree[i][j] = exp(-pow(10, 9) * pow(Y_tree[i][j] - delay, 2));
		}

	return J_tree;
}