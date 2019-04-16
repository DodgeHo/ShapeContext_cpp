#pragma once
# define PI 3.14159265358979323846
# define r_inner  0.125
# define r_outer  2
# define nbins_r  5
# define nbins_theta 12
#include <math.h>
#include <numeric>
#include <algorithm>
#include <armadillo>
using namespace arma;
using namespace std;

pair<mat, double> sc_compute(mat Bsamp, mat Tsamp, vector<bool> out_vec, double mean_dist=0.0);

mat dist2(mat x, mat c);
mat Atan2(mat Y, mat X);
mat rem(mat X, double y);

pair<mat, mat> bookstein(mat X, mat Y, double beta_k);