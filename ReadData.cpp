#include "stdafx.h"
#include "Neo_K_Means.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
void Neo_K_Means::ReadData(char* address) {
	size_type nEntry = N * dim;
	//Import data
	ifstream in(address);
	string c;
	char* stopstring;
	int pos = 0;
	
	cout.precision(15);
	//Read data.txt to make a data matrix
	for (size_type i = 0; i<nEntry; i++) {
		in >> c;
		data.mat[pos++] = strtold(c.c_str(), &stopstring);
	}
}
