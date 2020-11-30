#include "helper.h"

#include <vector>
using namespace std;

double sum_sub(vector<double> list, int start, int end) {
	double sum = 0;
	for (int i = start; i < end; i++) sum += list[i];
	return sum;
}