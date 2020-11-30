#include "rand.h"
#include <math.h>
#include <time.h>

using namespace std;


static   mt19937   mersenneTwister(time(0));

/*generate one uniform random variable (from 0 to 1)*/
double randUniform() {
    double ret = (mersenneTwister() + 0.5) / (mersenneTwister.max() + 1.0);
    return ret;
}

/*generate one bernoulli random variable (0 or 1)*/
double randBM() {
    double ret = (mersenneTwister() + 0.5) / (mersenneTwister.max() + 1.0);
    if (ret < 0.5) return -1;
    else return 1;
}

/*generate n bernoulli random variable (0 or 1)*/

vector<double> randBM(int n)
{
    vector<double> ret(n, 0.);
    for (int i = 0; i < n; i++)
        ret[i] = randBM();
    return ret;
}

/*generate one bernoulli random variable (exp(-lambda/Nstep)-1 or exp(-lambda/Nstep))*/
double randPoisson(double lambda, int Nstep) {
    double ret = (mersenneTwister() + 0.5) / (mersenneTwister.max() + 1.0);
    if (ret < exp(-lambda / Nstep)) return exp(-lambda / Nstep) - 1;
    else return exp(-lambda / Nstep);
}

/*generate n bernoulli random variable (exp(-lambda/Nstep)-1 or exp(-lambda/Nstep))*/
vector<double> randPoisson(int n, double lambda, int Nstep) {
    vector<double> ret(n, 0.);
    for (int i = 0; i < n; i++)
        ret[i] = randPoisson(lambda, Nstep);
    return ret;
}
