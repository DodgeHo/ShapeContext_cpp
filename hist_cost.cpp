#include "hist_cost.h"

// Function to calculate the histogram cost between two matrices
vector< vector<double> > hist_cost(mat BH1, mat BH2) {
	// Get the number of samples and bins from the matrices
	int nsamp = BH1.n_rows;
	int nbins = BH1.n_cols;

	// Normalize the matrices
	mat BH1n = BH1 / repmat(sum(BH1, 1) + epsn, 1, nbins);
	mat BH2n = BH2 / repmat(sum(BH2, 1) + epsn, 1, nbins);

	vector< vector<double> > HC(nsamp, vector<double>(nsamp));

	// Iterate through the samples and bins to calculate the histogram cost
	for (int i = 0; i < nsamp; ++i) {
		for (int j = 0; j < nsamp; ++j) {
			double sum = 0;
			for (int k = 0; k < nbins; ++k) {
				// Calculate the histogram cost for each bin
				sum += 0.5 * (BH1n(i, k) - BH2n(j, k)) * (BH1n(i, k) - BH2n(j, k))
					/ (BH1n(i, k) + BH2n(j, k) + epsn);
			}
			// Store the histogram cost in the vector of vectors
			HC[i][j] = sum;
		}
	}

	return HC;
}