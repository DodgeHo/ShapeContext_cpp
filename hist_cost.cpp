
#include "hist_cost.h"
vector< vector<double> > hist_cost(mat BH1, mat BH2) {
	int nsamp, nbins;
	nsamp = BH1.n_rows;
	nbins = BH1.n_cols;

	mat BH1n = BH1 / repmat(sum(BH1, 1) + epsn, 1, nbins);
	mat BH2n = BH2 / repmat(sum(BH2, 1) + epsn, 1, nbins);
	//mat HC=zeros(nbins,nbins);
	vector< vector<double> > HC(nsamp,vector<double>(nsamp));
	for (int i=0;i< nsamp;++i)
		for (int j = 0; j < nsamp; ++j)
			HC[i][j]=0;

	for (int i = 0;i < nsamp; ++i)
		for (int j = 0; j < nsamp; ++j) 
			for (int k = 0; k < nbins; ++k)
				HC[i][j]+= 0.5*(BH1n(i, k) - BH2n(j, k))* (BH1n(i, k) - BH2n(j, k))
								/ (BH1n(i, k) + BH2n(j, k) + epsn);
	
	return HC;
}

