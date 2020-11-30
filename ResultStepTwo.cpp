#include "ResultStepTwo.h"
#include "ResultStepOne.h"
#include "seasonality.h"
#include "helper.h"
#include <math.h>
#include <vector>
#include "rand.h"
using namespace std;



vector<double> P_StepTwo(int Nstep, int Nitera, int A, int K, double p0, double p1, double pi, int f0, int f1, double h1, double h2, double S_init, double bar_Alpha, double delay, double lambda, vector<double> pathCom, vector<double> pathInd, vector<double> pathJump, double dt, double sigmaCom, double sigmaInd, double*** phi_hat_tree, double*** psi_hat_tree, vector<double> season) {
	vector<double> Q_simulation(Nstep + 1, 0.);
	Q_simulation = Q(Nstep, pathCom, pathInd, dt, sigmaCom, sigmaInd, season);

	vector<double> Q_r_simulation(Nstep + 1, 0.);
	Q_r_simulation = Q_r(Nstep, pathCom, dt, sigmaCom, season);

	vector<double> bar_Q_simulation(Nstep + 1, 0.);
	bar_Q_simulation = bar_Q(Nstep, pathCom, dt, sigmaCom, season);

	vector<double> Y_simulation(Nstep + 1, 0.);
	Y_simulation = Y(Nstep, delay, pathJump, dt, lambda);

	vector<double> J_simulation(Nstep + 1, 0.);
	J_simulation = J(Nstep, Y_simulation, delay);

	vector<double> phi_hat_simulation(Nstep + 1, 0.);
	phi_hat_simulation = phi_hat(Nstep, phi_hat_tree, lambda, h2, Nitera, pathJump);

	vector<double> psi_hat_simulation(Nstep + 1, 0.);
	psi_hat_simulation = psi_hat(Nstep, psi_hat_tree, lambda, h1, Nitera, pathJump, pathCom);

	vector<double> bar_P_simulation(Nstep + 1, 0.);
	bar_P_simulation = bar_P(Nstep, p0, p1, K, pi, bar_Q_simulation, Q_r_simulation, f0, f1, bar_Alpha, dt, J_simulation, season);

	vector<double> S_hat_simulation(Nstep + 1, 0.);
	S_hat_simulation = S_hat(Nstep, phi_hat_simulation, psi_hat_simulation, bar_P_simulation, J_simulation, A, p1, f1, dt, S_init, pi, K);

	vector<double> alpha_hat_simulation(Nstep + 1, 0.);
	alpha_hat_simulation = alpha_hat(Nstep, phi_hat_simulation, S_hat_simulation, psi_hat_simulation, bar_P_simulation, J_simulation, A, p1, f1, pi, K);

	vector<double> P_simulation(Nstep + 1, 0.);
	for (int i = 0; i < Nstep + 1; i++) {
		P_simulation[i] = K * Q_simulation[i] + p0 + pi * p1 * Q_r_simulation[i] + (1 - pi) * p1 * (bar_Q_simulation[i] + alpha_hat_simulation[i]) + J_simulation[i] * (f0 + f1 * (bar_Q_simulation[i] - season[i] + alpha_hat_simulation[i] - bar_Alpha));
	}
	return P_simulation;

}

vector<double> phi(int Nstep, int h2, int C, int A, double dt) {
	vector<double> phi_simulation(Nstep + 1, 0.);
	phi_simulation[Nstep] = h2;
	for (int i = Nstep - 1; i > -1; i--) {
		phi_simulation[i] = phi_simulation[i + 1] + (C - pow(phi_simulation[i + 1], 2) / A) * dt;
	}
	return phi_simulation;
}

double psi(int timeK, int Nstep, int Nitera, int Nmc, int A, int K, double p0, double p1, double pi, int f0, int f1, double h1, double h2, double S_init, double bar_Alpha, double delay, double lambda, vector<double> pathCom, vector<double> pathInd, vector<double> pathJump, double dt, double sigmaCom, double sigmaInd, double*** phi_hat_tree, double*** psi_hat_tree, vector<double> phi_simulation, vector<double> season) {
	vector<double> P_simulation(Nstep + 1, 0.);
	vector<double> P_simulation1(Nstep + 1, 0.);
	vector<double> P_simulation2(Nstep + 1, 0.);
	vector<double> P_simulation3(Nstep + 1, 0.);
	double total = 0;
	double integralValue;
	double integralValue1;
	double integralValue2;
	double integralValue3;
	int currentTime;

	vector<double> pathComNew(Nstep, 0.);
	vector<double> pathIndNew(Nstep, 0.);
	vector<double> pathJumpNew(Nstep, 0.);


	for (int i = 0; i < Nmc; i++) {
		vector<double> pathComTail = randBM(Nstep - timeK);
		vector<double> pathIndTail = randBM(Nstep - timeK);
		vector<double> pathJumpTail = randPoisson(Nstep - timeK, lambda, Nstep);

		for (int j = 0; j < timeK; j++) {
			pathComNew[j] = pathCom[j];
			pathIndNew[j] = pathInd[j];
			pathJumpNew[j] = pathJump[j];
		}

		for (int j = 0; j < Nstep - timeK; j++) {
			pathComNew[j + timeK] = pathComTail[j];
			pathIndNew[j + timeK] = pathIndTail[j];
			pathJumpNew[j + timeK] = pathJumpTail[j];
		}

		P_simulation = P_StepTwo(Nstep, Nitera, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, pathComNew, pathIndNew, pathJumpNew, dt, sigmaCom, sigmaInd, phi_hat_tree, psi_hat_tree, season);

		integralValue = 0;
		for (int j = timeK; j < Nstep; j++) {
			integralValue += phi_simulation[j] / (A + K) * exp(-sum_sub(phi_simulation, timeK, j + 1) / (A + K) * dt) * P_simulation[j] * dt;
		}

		for (int j = 0; j < Nstep - timeK; j++) {
			pathComNew[j + timeK] = -pathComTail[j];
		}

		P_simulation1 = P_StepTwo(Nstep, Nitera, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, pathComNew, pathIndNew, pathJumpNew, dt, sigmaCom, sigmaInd, phi_hat_tree, psi_hat_tree, season);

		integralValue1 = 0;
		for (int j = timeK; j < Nstep; j++) {
			integralValue1 += phi_simulation[j] / (A + K) * exp(-sum_sub(phi_simulation, timeK, j + 1) / (A + K) * dt) * P_simulation1[j] * dt;
		}

		for (int j = 0; j < Nstep - timeK; j++) {
			pathComNew[j + timeK] = -pathComTail[j];
			pathIndNew[j + timeK] = -pathIndTail[j];
		}

		P_simulation2 = P_StepTwo(Nstep, Nitera, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, pathComNew, pathIndNew, pathJumpNew, dt, sigmaCom, sigmaInd, phi_hat_tree, psi_hat_tree, season);

		integralValue2 = 0;
		for (int j = timeK; j < Nstep; j++) {
			integralValue2 += phi_simulation[j] / (A + K) * exp(-sum_sub(phi_simulation, timeK, j + 1) / (A + K) * dt) * P_simulation2[j] * dt;
		}

		for (int j = 0; j < Nstep - timeK; j++) {
			pathComNew[j + timeK] = -pathComTail[j];
		}

		P_simulation3 = P_StepTwo(Nstep, Nitera, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, pathComNew, pathIndNew, pathJumpNew, dt, sigmaCom, sigmaInd, phi_hat_tree, psi_hat_tree, season);

		integralValue3 = 0;
		for (int j = timeK; j < Nstep; j++) {
			integralValue3 += phi_simulation[j] / (A + K) * exp(-sum_sub(phi_simulation, timeK, j + 1) / (A + K) * dt) * P_simulation3[j] * dt;
		}


		total += 0.25 * integralValue + 0.25 * integralValue1 + 0.25 * integralValue2 + 0.25 * integralValue3;
	}

	double mean = total / Nmc;
	return h1 * exp(-sum_sub(phi_simulation, timeK, Nstep + 1) / (A + K) * dt) - mean;
}

vector<double> S(int Nstep, vector<double> phi_simulation, vector<double> psi_simulation, vector<double> P_simulation, int A, int K, double dt, double S_init) {
	vector<double> AStar_simulation(Nstep + 1, 0.);
	for (int i = 0; i < Nstep + 1; i++) {
		AStar_simulation[i] = -(psi_simulation[i] + P_simulation[i]) / (A + K);
	}
	vector<double> integral1(Nstep + 1, 0.);
	for (int i = 1; i < Nstep + 1; i++) {
		integral1[i] = integral1[i - 1] + phi_simulation[i] / (A + K) * dt;
	}
	vector<double> integral2(Nstep + 1, 0.);
	for (int i = 1; i < Nstep + 1; i++) {
		integral2[i] = integral2[i - 1] + AStar_simulation[i] * exp(integral1[i]) * dt;
	}
	vector<double> S_simulation(Nstep + 1, 0.);
	S_simulation[0] = S_init;
	for (int i = 1; i < Nstep + 1; i++) {
		S_simulation[i] = exp(-integral1[i]) * (S_init + integral2[i]);
	}
	return S_simulation;
}


vector<double> alpha(int Nstep, vector<double> phi_simulation, vector<double> psi_simulation, vector<double> S_simulation, vector<double> P_simulation, int A, int K) {
	vector<double> alpha_simulation(Nstep + 1, 0.);
	for (int i = 0; i < Nstep + 1; i++) {
		alpha_simulation[i] = -(psi_simulation[i] + phi_simulation[i] * S_simulation[i] + P_simulation[i]) / (A + K);
	}
	return alpha_simulation;
}