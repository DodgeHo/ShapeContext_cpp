#include "EMG.h"


mat EMG(mat output) {

	int data_length, n;
	data_length = output.n_rows;
	n = output.n_cols;

	mat emgmat(n, n);
	emgmat.zeros();

	mat row_mat(data_length, 1);
	for (int i = 0; i < data_length; ++i)
		row_mat(i, 0) = i + 1.0;
	mat row_mat_join = join_rows(row_mat, output);

	HungarianAlgorithm HungAlgo;
	vector<int> cvec,cvec2;

	int k,s;
	double mse1, mean_dist,beta_k,cost;

	for (int i = 1; i <= n; ++i)
		for (int j = 0; j < i; ++j){
			mat X = join_rows(row_mat, output.col(i));
			mat Y = join_rows(row_mat, output.col(j));
	
			mat xk = X; //initialize transformed version of model pointset

			k = 1;//initialize counter
			s = 1;

			vector<bool> out_vec_1(data_length, false);
			vector<bool> out_vec_2(data_length, false);
			//out_vec_{1,2} are indicator vectors for keeping track of 
			//estimated outliers on each iteration

			while (s) {
				cout << "iter=" << k << endl;

				//compute shape contexts for (transformed) model
				pair< mat, double > SC1 = sc_compute(xk.t(), zeros(1, data_length), out_vec_1);
				mean_dist = SC1.second;

				//compute shape contexts for target, using the scale estimate from the warped model
				//Note: this is necessary only because out_vec_2 can change on each iteration,
				//which affects the shape contexts.
				//Otherwise,Y does not change.
				pair< mat, double > SC2 = sc_compute(Y.t(), zeros(1, data_length), out_vec_2, mean_dist);
								
				//compute regularization parameter
				beta_k = (mean_dist * mean_dist) * beta_init * pow( annealing_rate, (k - 1)) ;

				//compute pairwise cost between all shape contexts
				vector< vector<double> > costmat = hist_cost(SC1.first, SC2.first);

				/* when ndum>0
				if (ndum>0) {
					vector<double> temp(data_length);
					for (int j = 0; j < data_length; ++j){
						temp[j]=eps_dum;
					}

					for (int i = 0; i < ndum; ++i){
						costmat.push_back(temp);				
					}

					for (int j = 0; j < data_length+eps_dum; ++j){
						for (int i = 0; i < ndum; ++i){
							costmat[j].push_back(eps_dum);				
						}
					}
				}
				*/
				cvec2.clear();
				cout << "running hungarian alg." << endl;
				cost = HungAlgo.Solve(costmat, cvec2);
				costmat.clear();
				cout << "Hungarian alg done." << endl;
				//update outlier indicator vectors
				cvec.clear();
				cvec=BubbleSort(cvec2, data_length);
				for (int i = 0; i < data_length; i++)
					if (cvec2[i] > data_length)
						out_vec_1[i] = 1;
				for (int i = 0; i < data_length; i++)
					if (cvec[i] > data_length)
						out_vec_2[i] = 1;
				
				//ormat versions of Xk and Y that can be plotted with outliers'correspondences missing
				mat X2b = zeros(data_length,2);
				for (int i = 0; i < data_length; i++)
					X2b.row(i) = X.row(cvec[i]);
				

				//estimate regularized TPS transformation
				pair<mat, mat> cxy=bookstein(X2b, Y, beta_k);
				mat cx = cxy.first;
				mat cy = cxy.second;
		
				/*
				//calculate affine cost
				mat A = join_rows(cx.rows(data_length + 1,data_length + 2), 
								cy.rows(data_length + 1, data_length + 2));
				mat U, V;
				colvec s;
				svd(U, s, V, A);*/

				//warp each coordinate
				mat U = dist2(X2b, X);
				for (int i = 0; i < data_length; i++)
					for (int j = 0; j < data_length; j++)
						if (U(i, j) < 0)
							U(i, j) = 0;
				
				U = U % log(U + epsn);
			
				mat fx_aff = (cx.rows(data_length , data_length + 2).t())
					* join_cols(ones(1,data_length), X.t());

				
				mat fx_wrp = (cx.rows(0,data_length-1).t())*U;
				
				mat fy_aff = (cy.rows(data_length, data_length + 2).t())
					* join_cols(ones(1, data_length), X.t());
				mat fy_wrp = (cy.rows(0, data_length - 1).t()) * U;
				
				mat Z = join_cols(fx_aff + fx_wrp, fy_aff + fy_wrp).t();
			
				/*
				% compute the mean squared error between synthetic warped image
				% and estimated warped image (using ground-truth correspondences
				% on TPS transformed image)
				%mse2 = mean((Y(:,1)-Z(:,1)).^2 + (Y(:,2)-Z(:,2)).^2);
            
				% Chui actually does mean of non-squared distance
				*/
				U.clear();
				U = sqrt(
					( Y.col(0) - Z.col(0))% (Y.col(0) - Z.col(0)) 
					+ 
					(Y.col(1) - Z.col(1))% (Y.col(1) - Z.col(1))
					);
				
				mse1 = accumulate(U.begin(), U.end(), 0.0) * 1.0 / U.size();
				cout << "distance = " << mse1;



				// update Xk for the next iteration
					xk = Z;

				// stop early if shape context score is sufficiently low
					if (k == emg_iter_time)
						s = 0;
					else
						k++;
			}



			emgmat(i, j) = mse1;
		}


	for (int i = 0; i <= n; ++i)
		for (int j = i + 1; j <= n; ++j)
			emgmat(i, j) = emgmat(j, i);
	return emgmat;
}



vector<int> BubbleSort(vector<int> p, int length)
{
	vector<int> ind_diff(length);
	for (int m = 0; m < length; m++)
	{
		ind_diff[m] = m;
	}

	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length - i - 1; j++)
		{
			if (p[j] > p[j + 1])
			{
				int temp = p[j];
				p[j] = p[j + 1];
				p[j + 1] = temp;

				int ind_temp = ind_diff[j];
				ind_diff[j] = ind_diff[j + 1];
				ind_diff[j + 1] = ind_temp;
			}
		}
	}
	return ind_diff;
}

