#pragma once

#include <vector>
#include <random>

/*generate one uniform random variable (from 0 to 1)*/
double randUniform();

/*generate one bernoulli random variable (0 or 1)*/
double randBM();

/*generate n bernoulli random variable (0 or 1)*/
std::vector<double> randBM(int n);

/*generate one bernoulli random variable (exp(-lambda/Nstep)-1 or exp(-lambda/Nstep))*/
double randPoisson(double lambda, int Nstep);

/*generate n bernoulli random variable (exp(-lambda/Nstep)-1 or exp(-lambda/Nstep))*/
std::vector<double> randPoisson(int n, double lambda, int Nstep);