#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include "Neo_K_Means.h"

using namespace std;

void Neo_K_Means::Estimate_alpha_beta(double alpha_delta, double beta_delta) {

	size_type clusterNumber, alphaN, betaN, count;
	Matrix<double> D(N, K), C(K, dim), temp(1, dim);
	vector<pair<double, size_type> > pv(N);
	vector<double> total_sum_by_cluster(K), min_sum_by_cluster(K), square_sum_by_cluster(K);
	vector<int> count_cluster(K);
	double min, mean, std, total_sum = 0, min_sum = 0, square_sum = 0, criterion;
	//Calculate the coordinations of mean and save them in a matrix 'C'.
	for (size_type j = 0; j < K; j++) {
		count = 0;
		for (size_type i = 0; i < N; i++) {
			if (initU(i, j)) {
				temp += data(i, ':');
				count++;				
			}
		}
		for (size_type i = 0; i < dim; i++) {
			C(j,i) = temp((size_type)0, i) / count;
		}

		temp.mat.assign(dim, 0);
	}
	C.ShowInfo(1);
	for (size_type i = 0; i<N; i++) {
		min = numeric_limits<double>::infinity();
		for (size_type j = 0; j<K; j++) {
			if (min >(D(i, j) = sqrt((data(i, ':') - C(j, ':')).Norm()))) {
				min = D(i, j);
				clusterNumber = j;
			}
			total_sum_by_cluster[j] += D(i, j);
			total_sum += D(i, j);
		}
		pv[i] = (pair<double, size_type>(min, clusterNumber));

		min_sum += min;
		square_sum += min *min;
		min_sum_by_cluster[clusterNumber] += min;
		square_sum_by_cluster[clusterNumber] += min*min;
		count_cluster[clusterNumber]++;
	}


	//Estimate beta
	betaN = 0;
	mean = min_sum / N;
	std = sqrt(square_sum / N - mean*mean);
	criterion = mean + beta_delta*std;
	for (size_type i = 0; i<N; i++)
		if (pv[i].first > criterion) betaN++;
	beta = (double)betaN / N;

	//Estimate alpha
	alphaN = 0;
	if (alpha_delta == 0) {
		criterion = (1 / (1 + K));
		for (size_type i = 0; i<K; i++) {
			for (size_type j = 0; j < N; j++)
				if ((D(j, i) / total_sum_by_cluster[i]) < criterion) alphaN++;
		}
		alpha = (double)alphaN / N;
	}

	else {
		bool flag = false;
		for (size_type j = 0; j<K; j++) {
			for (size_type i = 0; i<K; i++) {
				if (count_cluster[i] != 0) {
					flag = true;
					break;
				}
			}
			if (!flag) continue;

			mean = min_sum_by_cluster[j] / count_cluster[j];
			std = sqrt(square_sum_by_cluster[j] / count_cluster[j] - mean*mean);
			criterion = mean + alpha_delta*std;

			for (size_type i = 0; i<N; i++)
				if (pv[i].second != j && D(i, j) <= criterion) alphaN++;
		}
		alpha = (double)alphaN / N;

	}
	return;
}
