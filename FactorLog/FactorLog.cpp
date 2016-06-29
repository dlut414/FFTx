// FactorLog.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <bitset>
#include "../common/Complex.h"
using namespace std;

#define BIT 16

typedef double T;
typedef Complex<T> C;

template <typename R>
__forceinline C EI(const R& k) {
	return C(cos(k), sin(k));
}

void FactorLog(unsigned nn) {
	vector<C> factors;
	ofstream file;
	file.open("factors.txt");
	const T PI = T(3.14159265358979323846);
	for (unsigned step = 1; step < (nn); step <<= 3) {
		for (unsigned factorId = 0; factorId < step; factorId++) {
			const T theta = -(PI / T(step << 2))* factorId;
			const C factor1 = EI(theta);
			const C factor2 = EI(theta * 2);
			const C factor3 = EI(theta * 3);
			const C factor4 = EI(theta * 4);
			const C factor5 = EI(theta * 5);
			const C factor6 = EI(theta * 6);
			const C factor7 = EI(theta * 7);
			factors.push_back(factor1);
			factors.push_back(factor2);
			factors.push_back(factor3);
			factors.push_back(factor4);
			factors.push_back(factor5);
			factors.push_back(factor6);
			factors.push_back(factor7);
		}
	}
	for (unsigned i = 0; i < factors.size(); i++) {
		short imag, real;
		imag = short(factors[i].imag * 0x7fff);
		real = short(factors[i].real * 0x7fff);
		file << setfill('0') << setw(4) << hex << imag << setfill('0') << setw(4) << hex << real << "," << endl;
	}
	cout << factors.size() << std::endl;
	file.close();
}

int main()
{
	FactorLog(2 << 15);
    return 0;
}

