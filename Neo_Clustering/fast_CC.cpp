#include "stdafx.h"
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include "Neo_K_Means.h"

using namespace std;

size_t cal_count;

void Neo_K_Means::Calculate_center(Matrix<double> data, Matrix<double>& Center, vector<int> CenterNumber, vector<int> OtherCenterNumber, int version) {
	size_t A, B, C, D;
	//Calculate the distance of row
	if (version == 0) {
		A = N;
		B = dim;
		C = K;
		D = L;
	}
	//Calculate the distance of col
	else if (version == 1) {
		A = dim;
		B = N;
		C = L;
		D = K;
	}
	Matrix<double> temp_sum(1, B);
	int count1, count2;
	double _sum, avg;
	for (size_t c = 0; c < C; c++) {
		count1 = 0;
		fill(temp_sum.mat.begin(), temp_sum.mat.end(), 0);
		for (size_t a = 0; a < A; a++) {
			if (CenterNumber[a] == c) {
				temp_sum += data(a, ':');
				count1++;
			}
		}
		for (size_t d = 0; d < D; d++) {
			count2 = 0;
			_sum = 0;
			for (size_t b = 0; b < B; b++) {
				if (OtherCenterNumber[b] == d) {
					_sum += temp_sum((size_t)0, b);
					count2++;
				}
				else continue;
			}
			if (count2 == 0) continue;
			avg = _sum / (count1*count2);

			for (size_t b = 0; b < B; b++) {
				if (OtherCenterNumber[b] == d) {
					Center(c, b) = avg;
				}
			}
		}
	}
};

//Calculate distance between old center and current center and distance among current centers.
//Return minimu distance.
double CalculateDistanceOfCenter(Matrix<double> oldCenter, Matrix<double> Center, vector<double>& oldC2C, Matrix<double>& C2C) {
	size_t K = Center.Get_r();
	size_t dim = Center.Get_r();
	double min = numeric_limits<double>::infinity();

	for (size_t i = 0; i < K; i++) {
		oldC2C[i] = sqrt((Center(i, ':') - oldCenter(i, ':')).Norm());
	}

	for (size_t i = 0; i < K; i++) {
		for (size_t j = 0; j < K; j++) {
			if (i == j) continue;
			C2C(i, j) = sqrt((Center(i, ':') - Center(j, ':')).Norm());
			C2C(j, i) = C2C(i, j);
			if (min > C2C(i, j)) {
				min = C2C(i, j);
			}
		}
	}
	return min;
}

void UpdateBound(std::vector<double>& Upper, Matrix<double>& Lower, std::vector<int> CenterNumber, std::vector<double> oldC2C) {
	size_t N = Lower.Get_r(), K = Lower.Get_c();
	for (size_t i = 0; i < N; i++) {
		int temp = CenterNumber[i];
		Upper[i] = Upper[i] + oldC2C[temp];
		for (size_t j = 0; j < K; j++) {
			Lower(i, j) = Lower(i, j) - oldC2C[j]> 0 ? Lower(i, j) - oldC2C[j] : 0;
		}
	}
}


void Neo_K_Means::fast_CC(int iteration, double threshold, Matrix<int> initU, Matrix<int> initV) {

	size_t t = 0;
	/* size_t alphaN_row = (size_t)floor((alpha_row*N) + 0.5);
	size_t alphaN_col = (size_t)floor((alpha_col*dim) + 0.5);
	size_t betaN_row = floor((beta_row*N) + 0.5);
	size_t betaN_col = floor((beta_col*dim) + 0.5);
	*/
	double J = numeric_limits<double>::infinity();
	double oldJ = 0;

	vector<bool> flag_row(N), flag_col(dim);
	vector<int> CenterNumberOfRow(N), CenterNumberOfCol(dim);
	vector<double> J_track, oldC2C_row(K), oldC2C_col(L), Upper_row(N), Upper_col(dim);
	Matrix<double> C2C_row(K, K), C2C_col(L, L), D_row(N, K), D_col(dim, L);
	Matrix<double> Lower_row(N, K), Lower_col(dim, L);
	Matrix<double> Center_row(K, dim), Center_col(L, N), oldCenter_row(N, K), oldCenter_col(dim, L);
	double min, temp;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < K; j++) {
			if (initU(i, j) == 1) CenterNumberOfRow[i] = j;
		}
	}

	for (size_t i = 0; i < dim; i++) {
		for (size_t j = 0; j < L; j++) {
			if (initV(i, j) == 1) CenterNumberOfCol[i] = j;
		}
	}
	Calculate_center(data, Center_row, CenterNumberOfRow, CenterNumberOfCol, 0);
	Calculate_center(dataT, Center_col, CenterNumberOfCol, CenterNumberOfRow, 1);


	cout << "********** NEO-CC **********" << endl;
	while (abs(J - oldJ) > threshold && t <= iteration) {
		cal_count = 0;
		oldJ = J;
		J = 0;
		oldCenter_row = Center_row;
		oldCenter_col = Center_col;
		//--------------------------update row cluster--------------------------------------
		//Calculate the position of centers
		Calculate_center(data, Center_row, CenterNumberOfRow, CenterNumberOfCol, 0);

		//Calculate the mean
		min = CalculateDistanceOfCenter(oldCenter_row, Center_row, oldC2C_row, C2C_row);


		//Update upper and lower bound
		UpdateBound(Upper_row, Lower_row, CenterNumberOfRow, oldC2C_row);

		//Calculate the distance
		min = 1 / 2 * min;
		for (size_t i = 0; i < N; i++) {
			temp = sqrt((data(i, ':') - Center_row((size_t)CenterNumberOfRow[i], ':')).Norm());
			D_row(i, (size_t)CenterNumberOfRow[i]) = temp;
			Lower_row(i, (size_t)CenterNumberOfRow[i]) = temp;
			if (Upper_row[i] < min) {
				flag_row[i] = true;
			}
		}
		cal_count += N;
		for (size_t k = 0; k < K; k++) {
			for (size_t i = 0; i < N; i++) {
				if (flag_row[i] || CenterNumberOfRow[i] == k) continue;
				if (D_row(i, (size_t)CenterNumberOfRow[i]) > Lower_row(i, k) && D_row(i, (size_t)CenterNumberOfRow[i]) > 1 / 2 * C2C_row((size_t)CenterNumberOfRow[i], k)) {
					temp = sqrt((data(i, ':') - Center_row(k, ':')).Norm());
					Lower_row(i, k) = temp;
					D_row(i, k) = temp;
					cal_count++;
					if (temp < D_row(i, (size_t)CenterNumberOfRow[i])) {
						CenterNumberOfRow[i] = k;
						Upper_row[i] = temp;
					}
				}
			}
		}

		//--------------------------update col cluster--------------------------------------
		//Calculate the position of centers
		Calculate_center(dataT, Center_col, CenterNumberOfCol, CenterNumberOfRow, 1);

		//Calculate the mean
		min = CalculateDistanceOfCenter(oldCenter_col, Center_col, oldC2C_col, C2C_col);

		//Update upper and lower bound
		UpdateBound(Upper_col, Lower_col, CenterNumberOfCol, oldC2C_col);


		//Calculate the distance
		min = 1 / 2 * min;
		for (size_t i = 0; i < dim; i++) {
			temp = sqrt((dataT(i, ':') - Center_col((size_t)CenterNumberOfCol[i], ':')).Norm());
			D_col(i, (size_t)CenterNumberOfCol[i]) = temp;
			Lower_col(i, (size_t)CenterNumberOfCol[i]) = temp;
			if (Upper_col[i] < min) {
				flag_col[i] = true;
			}
		}
		cal_count += dim;
		for (size_t k = 0; k < L; k++) {
			for (size_t i = 0; i < dim; i++) {
				if (flag_col[i] && CenterNumberOfCol[i] == k) continue;
				if (D_col(i, (size_t)CenterNumberOfCol[i]) > Lower_col(i, k) && D_col(i, (size_t)CenterNumberOfCol[i]) > 1 / 2 * C2C_col((size_t)CenterNumberOfCol[i], k)) {
					temp = sqrt((dataT(i, ':') - Center_col(k, ':')).Norm());
					Lower_col(i, k) = temp;
					D_col(i, k) = temp;
					cal_count++;
					if (temp < D_col(i, (size_t)CenterNumberOfCol[i])) {
						CenterNumberOfCol[i] = k;
						Upper_col[i] = temp;
					}
				}
			}
		}

		for (size_t i = 0; i < N; i++) {
			J += D_row(i, (size_t)CenterNumberOfRow[i]) * D_row(i, (size_t)CenterNumberOfRow[i]);
		}
		for (size_t i = 0; i < dim; i++) {
			J += D_col(i, (size_t)CenterNumberOfCol[i]) * D_col(i, (size_t)CenterNumberOfCol[i]);
		}

		//increase # of iteration && record the value of objective function.      
		t++;
		J_track.push_back(J);
		cout << "***** iteration : " << t << ", objective : " << J << endl;
		cout << "***** calculation_count : " << cal_count << endl;
	}

	cout << "***** No. of iterations done: " << t << endl;

	cout << "***** (ROW) Total no. of data points: " << N << endl;
	//cout << "***** (ROW) alpha_r: " << alpha_row << ", alphaN_r: " << alphaN_row << endl;
	//cout << "***** (ROW) beta_r: " << beta_row << ", betaN_r: " << betaN_row << endl;
	cout << "***** (COL) Total no. of data points: " << dim << endl;
	//cout << "***** (COL) alpha_c: " << alpha_col << ", alphaN_c: " << alphaN_col << endl;
	//cout << "***** (COL) beta_c: " << beta_col << ", betaN_c: " << betaN_col << endl;
	return;

}
