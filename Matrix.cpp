#include "stdafx.h"
#include <iostream>
#include <vector>
#include "Neo_K_Means.h"
using namespace std;

template<typename T>
Matrix<T>::Matrix() {};

template<class T>
Matrix<T>::Matrix(int _r, int _c) : r(_r), c(_c)
{
	mat.resize(r*c);
}

template<class T>
Matrix<T>::Matrix(std::vector<T>& arr, int _r, int _c) : r(_r), c(_c), mat(arr)
{
}

template<class T>
int Matrix<T>::Get_r() const { return r; }
template<class T>
int Matrix<T>::Get_c() const { return c; }

template<class T>
void Matrix<T>::ShowInfo(int option) {
	cout << "r, c ::" << r << ", " << c << endl;
	cout << "# of matrix entry ::" << mat.size() << endl;
	if (option == 1) {
		for (size_t i = 0; i<r; i++) {
			size_t temp = i*c;
			for (size_t j = 0; j<c; j++) {
				cout << mat[temp + j] << ' ';
			}
			cout << endl << endl;
		}
	}
}

template<class T>
Matrix<T>& Matrix<T>::operator +=(const Matrix<T>& M) {
	//Check whether sizes are the same.
	if (r != M.Get_r() || c != M.Get_c()) {
		cout << "The sizes of matrix doesn't match." << endl;
		return *this;
	}

	int temp;
	for (int i = 0; i<r; i++) {
		temp = i*c;
		for (int j = 0; j<c; j++) {
			mat[temp + j] += M.mat[temp + j];
		}
	}

	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator +(const Matrix<T>&M) {
	Matrix<T> temp = *this;
	return temp += M;
}

template<class T>
Matrix<T>& Matrix<T>::operator -=(const Matrix<T>& M) {
	//Check whether sizes are the same.
	if (r != M.Get_r() || c != M.Get_c()) {
		cout << "The sizes of matrix doesn't match." << endl;
		return *this;
	}

	int temp;
	for (int i = 0; i<r; i++) {
		temp = i*c;
		for (int j = 0; j<c; j++) {
			mat[temp + j] -= M.mat[temp + j];
		}
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator -(const Matrix<T>&M) {
	Matrix<T> temp = *this;
	return temp -= M;
}

template<class T>
Matrix<T>& Matrix<T>::operator =(const Matrix<T>&M) {
	r = M.Get_r();
	c = M.Get_c();
	mat = M.mat;
	return *this;
}

template<class T>
T& Matrix<T>::operator ()(const size_t a, const size_t b) {
	//check the range of an given index
	if (a >= r || b >= c) {
		cout << "The given index 'a,b' is out of the range <operator (a,b) function>." << endl;
	
	}
	return mat[a*c + b];
}

template<class T>
Matrix<T> Matrix<T>::operator ()(const char a, const size_t b) {
	//check the range of an given index
	if (b >= c) {
		cout << "The given index 'b' is out of the range <operator (a,b) function>." << endl;
	}

	Matrix temp(r, 1);
	temp.mat.resize(r);
	for (int i = 0; i<r; i++) {
		temp.mat[i] = mat[i*c + b];
	}
	return temp;
}

template<class T>
Matrix<T> Matrix<T>::operator ()(const size_t a, const char b) {
	//check the range of an given index
	if (a >= r) {
		cout << "The given index 'a' is out of the range <operator (a,b) function>." << endl;
	}
	Matrix<T> temp(1, c);
	temp.mat.assign(mat.begin() + a*c, mat.begin() + (a + 1)*c);
	return temp;
}

template<class T>
double Matrix<T>::Norm() {
	double sum = 0;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j<c; j++) {
			sum += (mat[i*c + j] * mat[i*c + j]);
		}
	}
	return sum;
}

template class Matrix<double>;
template class Matrix<int>;
