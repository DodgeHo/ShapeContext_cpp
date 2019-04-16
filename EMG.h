#pragma once

#define emg_iter_time 5
#define beta_init 1 //initial regularization parameter (normalized)
#define annealing_rate 1 //annealing rate
#define ndum 0
#define eps_dum  0.15
#include "sc_compute.h"
#include "hist_cost.h"
#include "hungarian.h"

using namespace arma;
mat EMG(mat output);
vector<int> BubbleSort(vector<int> p, int length);