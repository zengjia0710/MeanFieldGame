#include "seasonality.h"
#include "ResultStepOne.h"
#include "ResultStepTwo.h"
#include <vector>
#include "rand.h"
#include <math.h>
using namespace std;

double OptimalValue(int Nstep, int Nmc, int Nitera, vector<double> season, vector<double> phi_simulation, int A, int C, int K, double p0, double p1, double bar_Alpha, int f0, int f1, int h0, int h1, int h2, double dt, double pi, double S_init, double delay, double lambda, double sigmaCom, double sigmaInd, double*** phi_hat_tree, double*** psi_hat_tree) {
	double total = 0;
	double integralValue;

	vector<double> randnCom;
	vector<double> randnInd;
	vector<double> randnPoi;

	vector<double> Q_r_simulation(Nstep + 1, 0.);
	vector<double> bar_Q_simulation(Nstep + 1, 0.);
	vector<double> Q_simulation(Nstep + 1, 0.);
	vector<double> Y_simulation(Nstep + 1, 0.);
	vector<double> J_simulation(Nstep + 1, 0.);

	vector<double> phi_hat_simulation(Nstep + 1, 0.);
	vector<double> psi_hat_simulation(Nstep + 1, 0.);
	vector<double> bar_P_simulation(Nstep + 1, 0.);
	vector<double> S_hat_simulation(Nstep + 1, 0.);
	vector<double> alpha_hat_simulation(Nstep + 1, 0.);

	vector<double> price(Nstep + 1, 0.);

	vector<double> psi_simulation(Nstep + 1, 0.);
	psi_simulation[Nstep] = h1;
	vector<double> P_simulation(Nstep + 1, 0.);
	vector<double> S_simulation(Nstep + 1, 0.);
	vector<double> alpha_simulation(Nstep + 1, 0.);

	for (int i = 0; i < Nmc; i++) {
		randnCom = randBM(Nstep);
		randnInd = randBM(Nstep);
		randnPoi = randPoisson(Nstep, lambda, Nstep);


		Q_r_simulation = Q_r(Nstep, randnCom, dt, sigmaCom, season);
		bar_Q_simulation = bar_Q(Nstep, randnCom, dt, sigmaCom, season);
		Q_simulation = Q(Nstep, randnCom, randnInd, dt, sigmaCom, sigmaInd, season);
		Y_simulation = Y(Nstep, delay, randnPoi, dt, lambda);
		J_simulation = J(Nstep, Y_simulation, delay);


		phi_hat_simulation = phi_hat(Nstep, phi_hat_tree, lambda, h2, Nitera, randnPoi);
		psi_hat_simulation = psi_hat(Nstep, psi_hat_tree, lambda, h1, Nitera, randnPoi, randnCom);
		bar_P_simulation = bar_P(Nstep, p0, p1, K, pi, bar_Q_simulation, Q_r_simulation, f0, f1, bar_Alpha, dt, J_simulation, season);
		S_hat_simulation = S_hat(Nstep, phi_hat_simulation, psi_hat_simulation, bar_P_simulation, J_simulation, A, p1, f1, dt, S_init, pi, K);
		alpha_hat_simulation = alpha_hat(Nstep, phi_hat_simulation, S_hat_simulation, psi_hat_simulation, bar_P_simulation, J_simulation, A, p1, f1, pi, K);


		for (int i = 0; i < Nstep + 1; i++) {
			price[i] = p0 + 87.4286117 * pi * Q_r_simulation[i] + 87.4286117 * (1 - pi) * (bar_Q_simulation[i] + alpha_hat_simulation[i]);
		}


		for (int j = 0; j < Nstep; j++) {
			psi_simulation[j] = psi(j, Nstep, Nitera, Nmc, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, randnCom, randnInd, randnPoi, dt, sigmaCom, sigmaInd, phi_hat_tree, psi_hat_tree, phi_simulation, season);
		}


		P_simulation = P_StepTwo(Nstep, Nitera, A, K, p0, p1, pi, f0, f1, h1, h2, S_init, bar_Alpha, delay, lambda, randnCom, randnInd, randnPoi, dt, sigmaCom, sigmaInd, phi_hat_tree, psi_hat_tree, season);
		S_simulation = S(Nstep, phi_simulation, psi_simulation, P_simulation, A, K, dt, S_init);
		alpha_simulation = alpha(Nstep, phi_simulation, psi_simulation, S_simulation, P_simulation, A, K);

		integralValue = 0;
		for (int k = 0; k < Nstep + 1; k++) {
			integralValue += (1 - pi) * (0.5 * A * pow(alpha_simulation[k], 2) + 0.5 * C * pow(S_simulation[k], 2)) * dt;
			integralValue += (1 - pi) * (0.5 * K * pow(Q_simulation[k] + alpha_simulation[k], 2)) * dt;
			integralValue += (1 - pi) * ((Q_simulation[k] + alpha_simulation[k]) * price[k]) * dt;
			integralValue += (1 - pi) * (J_simulation[k] * (Q_simulation[k] - bar_Q_simulation[k] + alpha_simulation[k] - bar_Alpha) * (f0 + 10000 * (bar_Q_simulation[k] - season[k] + alpha_hat_simulation[k] - bar_Alpha))) * dt;
			integralValue += pi * (Q_r_simulation[k] * price[k] + 0.5 * K * pow(Q_r_simulation[k], 2)) * dt;
		}
		total += integralValue + (1 - pi) * (h0 + h1 * S_simulation[Nstep] + h2 * pow(S_simulation[Nstep], 2));
	}
	double mean = total / Nmc;
	return mean;
}