// fftLab.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES
#include "FFTx.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

#define BIT 12
#define POINTS (POWER2<BIT>::value)
#define FREQ 50.0
typedef float TYPE;

DFT<TYPE, BIT> dft;
Radix_2<TYPE, BIT> radix2;
Radix_4<TYPE, BIT> radix4;
Radix_8<TYPE, BIT> radix8;

using namespace std;
int main() {
	const TYPE PI = TYPE(3.14159265358979323846);
	///test 1
	for (int i = 0; i < POINTS; i++) {
		dft.x[i] = Complex<TYPE>(sin(i*FREQ*2*PI / (POINTS-1)), 0);
	}
	///test 2
	//for (int i = 0; i < POINTS; i++) {
	//	if (i < POINTS >> 1) {
	//		dft.x[i] = Complex<TYPE>(1, 0);
	//	}
	//	else {
	//		dft.x[i] = Complex<TYPE>(0, 0);
	//	}
	//}
	for (int i = 0; i < POINTS; i++) {
		radix2.x[i] = dft.x[i];
		radix4.x[i] = dft.x[i];
		radix8.x[i] = dft.x[i];
	}

	dft.Run();
	radix2.Run_DIF_aligned();
	radix4.Run_DIF_aligned();
	radix8.Run_DIF_aligned();

	cout << " input " << endl;
	ofstream inputOut;
	inputOut.open("./inputOut.txt");
	//cout << " # \t Time \t Real \t Imag " << endl;
	for (int i = 0; i < POINTS; i++) {
		//cout << i << fixed << "\t" << i << "\t" << dft.x[i].real << "\t" << dft.x[i].imag << endl;
		inputOut << i << fixed << "\t" << i << "\t" << dft.x[i].real << "\t" << dft.x[i].imag << endl;
	}
	inputOut.close();

	cout << " output dft: " << endl;
	ofstream dftout;
	dftout.open("./dftout.txt");
	for (int i = 0; i < POINTS; i++) {
		//cout << i << fixed << "\t" << dft.out[i].real << "\t" << dft.out[i].imag << endl;
		dftout << i << fixed << "\t" << dft.out[i].real << "\t" << dft.out[i].imag << endl;
	}
	dftout.close();

	cout << " output radix-2: " << endl;
	ofstream rad2out;
	rad2out.open("./rad2out.txt");
	for (int i = 0; i < POINTS; i++) {
		//cout << i << fixed << "\t" << radix2.x[i].real << "\t" << radix2.x[i].imag << endl;
		rad2out << i << fixed << "\t" << radix2.x[i].real << "\t" << radix2.x[i].imag << endl;
	}
	rad2out.close();

	cout << " output radix-4: " << endl;
	ofstream rad4out;
	rad4out.open("./rad4out.txt");
	for (int i = 0; i < POINTS; i++) {
		//cout << i << fixed << "\t" << radix4.x[i].real << "\t" << radix4.x[i].imag << endl;
		rad4out << i << fixed << "\t" << radix4.x[i].real << "\t" << radix4.x[i].imag << endl;
	}
	rad4out.close();

	cout << " output radix-8: " << endl;
	ofstream rad8out;
	rad8out.open("./rad8out.txt");
	for (int i = 0; i < POINTS; i++) {
		//cout << i << fixed << "\t" << radix8.x[i].real << "\t" << radix8.x[i].imag << endl;
		rad8out << i << fixed << "\t" << radix8.x[i].real << "\t" << radix8.x[i].imag << endl;
	}
	rad8out.close();

    return 0;
}
