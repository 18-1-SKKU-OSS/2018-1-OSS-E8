#include "stdafx.h"
#include "Neo_K_Means.h"
#include <string>

#include <iostream>
//#include "mat.h"
//#include "matrix.h"

using namespace std;

int main(int argc, char** argv) {
/*
	//Import data
	ifstream in("C:\\Users\\YS\\Dropbox\\coding\\initU.txt");
	string c;
	char* stopstring;
	int pos = 0;
	vector<int> readU(2000);
	//Read data.txt to make a data matrix
	for (size_type i = 0; i<2000; i++) {
		in >> c;
		readU[pos++] = atoi(c.c_str());
	}
	for (size_type i = 0; i < 2000; i++) {
		if (readU[i] != algo.initU.mat[i]) {
			cout << i << " : different" << endl;
		}
	}
	*/

	//Call the Neo_K_Means
	Neo_K_Means algo(2, 1000, 2);
	algo.ReadData("C:\\Users\\YS\\Dropbox\\coding\\data.txt");
	algo.Initialize();
	algo.Estimate_alpha_beta(1.9, 6);
	algo.Neo_kmeans(1000, 0.00001);
	
	return 0;
}
