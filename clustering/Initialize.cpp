#include "stdafx.h"
#include <vector>
#include <random>
#include <limits>
#include "Neo_K_Means.h"

using namespace std;
void Neo_K_Means::Initialize() {
	//Select a set of seed using k-means++
	Matrix<int> initU(N, K);
	vector<double> min_distance_of_each_data(N);
	vector<int> cluster_of_each_data(N);
	vector<int> seedID(K);
	vector<bool> skipPoints(N, false);

	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(1, N);

	double temp, totalDistance = 0, ran;

	//choose first seed.
	seedID[0] = dis(gen) - 1;
	skipPoints[seedID[0]] = true;
	for (size_t i = 0; i < N; i++) {
		temp = (data(i, ':') - data(seedID[0], ':')).Norm();
		min_distance_of_each_data[i] = temp;
		totalDistance += temp;
	}

	for (size_t i = 1; i < K; i++) {
		//choose a new seed 
		uniform_real_distribution<> dis(0, totalDistance);
		ran = dis(gen);
		for (size_t j = 0; j < N; j++) {
			//Don't consider points already selected as starting points.
			if (!skipPoints[j]) {
				ran -= min_distance_of_each_data[j];
				if (ran <= 0) {
					//select this point as a new starting point
					seedID[i] = j;
					skipPoints[j] = true;
					break;
				}
			}
		}


		for (size_t j = 0; j < N; j++) {
			//calculate a distance between data and the new seed.
			temp = (data(j, ':') - data(seedID[i], ':')).Norm();

			//find smallest distance of each data from a new set of seed.
			if (min_distance_of_each_data[j] > temp) {
				totalDistance = totalDistance - min_distance_of_each_data[j] + temp;
				min_distance_of_each_data[j] = temp;
				cluster_of_each_data[j] = i;
			}
		}
	}

	for (size_t i = 0; i < N; i++) {
		initU(i, (size_t)cluster_of_each_data[i]) = 1;
	}
	this->initU = initU;

	//Do the traditional k-means algorithm by calling neo-k-means funtions with alpha = 0 and beta = 0.
	Neo_kmeans(1000, 0.001);
	this->initU = this->U;
	
}