// Neo-K-Means.cpp : 콘솔 응용 프로그램에 대한 진입점을 정의합니다.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <utility>
#include "Neo_K_Means.h"

using namespace std;

//This function is used when pair vector having elements of distance is sorted.
//Use decreasing order.
template <typename T>
bool comp(const pair<T, pair<size_type, size_type> >& a, const pair<T, pair<size_type, size_type> >& b)
{
	return a.first<b.first;
}

Neo_K_Means::Neo_K_Means(size_type _K, size_type _N, size_type _dim, double _alpha, double _beta)
	: K(_K), N(_N), dim(_dim), alpha(_alpha), beta(_beta), data(_N, _dim), U(_N, _dim)
{
}

void Neo_K_Means::Neo_kmeans(int iteration, double threshold) {

	Matrix<int> U = initU;

	size_type t = 0;
	size_type alphaN = (size_type)floor((alpha*N) + 0.5);
	size_type betaN = floor((beta*N) + 0.5);
	size_type numAssign = N - betaN;
	double J = numeric_limits<double>::infinity();
	double oldJ = 0;
	double epsilon = 0;
	vector<double> J_track;
	Matrix<double> D(N, K);

	while (abs(J - oldJ) > threshold && t <= iteration) {
		oldJ = J;
		J = 0;
		Matrix<double> mean(1, dim);
		//Compute cluster means and distances
		for (size_type j = 0; j<K; j++) {
			int sum = 0, count = 0;
			for (size_type i = 0; i<N; i++) {
				if (U(i, j) == 1) {
					mean += data(i, ':');
					count++;
				}
			}

			for (size_type i = 0; i<dim; i++) {
				mean(static_cast<size_type> (0), i) /= count;
			}

			for (size_type i = 0; i<N; i++) {
				D(i, j) = (data(i, ':') - mean).Norm();
			}
			mean.mat.assign(dim, 0);
		}

		//Sorting
		//vecotr<vector<T, size_type, size_type> > > pv, pv2;
		
		vector<pair<double, pair<size_type, size_type> > > pv, pv2;
		pv.reserve(N);			//pv is used for the first assignment of non-overlapped data.
		pv2.reserve(N*K);		//pv2 is used for the second assignment of data which is overlapped or not.

								//Make pv and pv2
		for (size_type i = 0; i<N; i++) {
			double min = numeric_limits<double>::infinity();
			size_type clusterNumber;
			for (size_type j = 0; j<K; j++) {
				pv2.push_back(pair<double, pair<size_type, size_type> >(D(i, j), pair<size_type, size_type>(i, j)));
				if (D(i, j) < min) {
					min = D(i, j);
					clusterNumber = j;
				}
			}
			pv.push_back(pair<double, pair<size_type, size_type> >(min, pair<size_type, size_type>(i, clusterNumber)));
		}

		sort(pv.begin(), pv.end(), comp<double>);


		// Make (N-betaN) assignments
		U.mat.assign(K*N, 0);
		size_type x, y;
		for (size_type i = 0; i<numAssign; i++) {
			J += pv[i].first;
			x = pv[i].second.first;
			y = pv[i].second.second;
			U(x, y) = 1;
			D(x, y) = numeric_limits<double>::infinity();
			pv2[x*K + y].first = numeric_limits<double>::infinity();
		}

		sort(pv2.begin(), pv2.end(), comp<double>);

		//Make (alphaN + betaN) assignments
		int count = 0;
		while (count < alphaN + betaN) {
			J += pv2[count].first;
			U(pv2[count].second.first, pv2[count].second.second) = 1;
			count++;
		}

		//increase # of iteration && record the value of objective function.		
		t++;
		J_track.push_back(J);
		cout << "***** iteration : " << t << ", objective : " << J << endl;
	}

	cout << "***** No. of iterations done: " << t << endl;
	this->U = U;	//The final indicator matrix is saved in member variable of neo_kmeans class. 

	cout << "***** Total no. of data points: " << N << endl;
	cout << "***** alpha: " << alpha << ", alphaN: " << alpha*N << endl;
	cout << "***** beta: " << beta << ", betaN: " << beta*N << endl;
	return;
}
