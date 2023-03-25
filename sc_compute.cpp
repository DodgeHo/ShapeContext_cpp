
#include "sc_compute.h"


pair<mat, double> sc_compute(mat Bsamp, mat Tsamp, vector<bool> out_vec, double mean_dist ) {
	/*		compute(r, theta) histograms for points along boundary
	Bsamp is 2 x data_length(xand y coords.)	Tsamp is 1 x data_length(tangent theta)
	out_vec is 1 x data_length(0 for inlier, 1 for outlier)

	mean_dist is the mean distance, used for length normalization
	if it is not supplied, then it is computed from the data

	outliers are not counted in the histograms, but they do get
	assigned a histogram
	*/
	int data_length;
	data_length = Bsamp.n_cols;
	//in_vec = (out_vec == 0);
	
	// compute r, theta arrays
	mat r_array = sqrt(dist2(Bsamp.t(), Bsamp.t()));

	// real is needed to prevent bug in Unix version
	mat theta_array = Atan2( Bsamp.row(1).t()*ones(1,data_length)-ones(data_length,1)*Bsamp.row(1),
	(Bsamp.row(0)).t() * ones(1, data_length) - ones(data_length, 1) * (Bsamp.row(0))).t();

	theta_array -= (Tsamp.t()*ones(1,data_length));

		// create joint(r, theta) histogram by binning r_array and
		// theta_array

		// normalize distance by mean, ignoring outliers
	if (!mean_dist) {
		/*
		vector<double> tmp;
		tmp.clear();
		for (int i = 0; i < out_vec.size(); ++i) {
			if (out_vec[i]==false)
				tmp.push_back(r_array(i, i));
		}	*/
		mat tmp = r_array;
		for (int i = 0; i < int(out_vec.size()); ++i) {
			if (out_vec[i]) {
				tmp.row(i) = 0;
				tmp.col(i) = 0;
			}
		}	
		//mean_dist = mean(rowvec(tmp));
		mean_dist = accumulate(tmp.begin(), tmp.end(),0.0)*1.0/tmp.size();
	}
	r_array /= mean_dist;

	//use a log. scale for binning the distances
	mat r_array_q = zeros(data_length, data_length);
	vec r_bin_edges = arma::logspace(log10(r_inner), log10(r_outer), nbins_r);
	for (auto r_value:r_bin_edges ){
		for (int i = 0; i < data_length; i++)
			for (int j = 0; j < data_length; j++)
				r_array_q (i,j)+= (r_array(i, j) < r_value);
	}

	
	// flag all points inside outer boundary
	mat fz = zeros(data_length, data_length);
	for (int i = 0; i < data_length; i++)
		for (int j = 0; j < data_length; j++)
			fz(i, j) = (r_array(i, j)<r_outer);
	
	// put all angles in[0, 2PI) range
	// quantize to a fixed set of angles(bin edges lie on 0, (2 * PI) / k, ...2 * PI
	mat	theta_array_q = rem(rem(theta_array, 2 * PI) + 2 * PI, 2 * PI);
	theta_array_q = floor(theta_array_q / (2 * PI / nbins_theta));

	int nbins = nbins_theta * nbins_r;
	int num;
	mat BH(data_length, nbins);
	BH.zeros();
	for (int n = 0; n < data_length; ++n) {
		mat fzn = fz.row(n);
		for (int i = 0; i <  int(out_vec.size()); ++i) {
			if (out_vec[i])
				fzn(0, i) = 0;
		}
		for (int i = 0; i < data_length; ++i) 
			if(fzn(0,i)){
				num = ((int)theta_array_q(n, i) + nbins_theta * (int)(r_array_q(n, i)-1));
				BH(n, num)++;
			}
	}
	
	return make_pair(BH, mean_dist);

}




mat dist2(mat x, mat c) {
	/*
	Description
	D = DIST2(X, C) takes two matrices of vectors and calculates the
	squared Euclidean distance between them.  Both matrices must be of
	the same column dimension.  If X has M rows and N columns, and C has
	L rows and N columns, then the result has M rows and L columns.  The
	I, Jth entry is the  squared distance from the Ith row of X to the
	Jth row of C.

	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
	*/
	int ndata, ncentres, dim;
	ndata = x.n_rows;
	ncentres = c.n_rows;
	dim = x.n_cols;
	mat A = (ones(ncentres, 1) * sum((x % x).t(), 0)).t()
		+ ones(ndata, 1) * sum((c % c).t(), 0)
		- (x * (c.t()) * 2);
	return A;
}

mat Atan2(mat Y, mat X) {
	int data_length, nbins;
	data_length = Y.n_rows;
	nbins = Y.n_cols;
	mat result(data_length, nbins);
	for (int i = 0; i < data_length; i++)
		for (int j = 0; j < nbins; j++)
			result(i, j) = atan2(Y(i, j), X(i, j));
	return result;
}

mat rem(mat X, double y) {
	int data_length, nbins;
	data_length = X.n_rows;
	nbins = X.n_cols;

	mat result(data_length, nbins);
	for (int i = 0; i < data_length; i++)
		for (int j = 0; j < nbins; j++)
			result(i, j) = X(i, j) - y * (int)(X(i, j) / y);
	return result;
}

pair<mat,mat> bookstein(mat X,mat Y,double beta_k)
{
	// Bookstein 
	int N = X.n_rows;
	// compute distances between left points
	mat r2 = dist2(X,X);
	mat K = r2%(log(r2+ eye(N, N))); // add identity matrix to make K zero on the diagonal
	mat P = join_rows(ones(N,1),X);
	mat L= join_rows(P.t(),zeros(3,3));

	P = join_rows(K, P);
	L = join_cols(P, L);
	//actually: L = [K  P
    //				 P' zeros(3,3)];
	mat V = join_rows(Y.t(),zeros(2,3));

	// regularization
	for (int i = 0; i < N; i++)
		L(i, i) += beta_k;

	mat c = (L.i()) * (V.t());
	
	mat cx = c.col(0);
	mat cy = c.col(1);
	return make_pair(cx, cy);
}