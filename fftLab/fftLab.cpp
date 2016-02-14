// fftLab.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES
#include "FFTx.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>

#define BIT 6
#define POINTS (POWER2<BIT>::value)
#define FREQ 1
#define DT (1./FREQ)
#define OMIGA (0.5*M_PI)
typedef float TYPE;

using namespace std;
int main()
{
	DFT<TYPE, BIT> dft;
	Radix_2<TYPE, BIT> radix2;
	Radix_4<TYPE, BIT> radix4;
	Radix_8<TYPE, BIT> radix8;
	for (int i = 0; i < POINTS; i++) {
		dft.x[i] = Complex<TYPE>(sin(OMIGA*i*DT), 0);
		radix2.x[i] = Complex<TYPE>(sin(OMIGA*i*DT), 0);
		radix4.x[i] = Complex<TYPE>(sin(OMIGA*i*DT), 0);
		radix8.x[i] = Complex<TYPE>(sin(OMIGA*i*DT), 0);
	}

	cout << " input " << endl;
	cout << " # \t Time \t Real \t Imag " << endl;
	for (int i = 0; i < POINTS; i++) {
		cout << i << fixed << setw(15) << i*DT << setw(15) << dft.x[i].real << " , " << setw(15) << dft.x[i].imag << endl;
	}

	dft.Run();
	radix2.Run_aligned();
	radix4.Run_aligned();
	radix8.Run_aligned();

	cout << " output dft: " << endl;
	for (int i = 0; i < POINTS; i++) {
		cout << i << fixed << setw(15) << dft.out[i].real << setw(15) << dft.out[i].imag << endl;
	}
	cout << " output radix-2: " << endl;
	for (int i = 0; i < POINTS; i++) {
		cout << i << fixed << setw(15) << radix2.x[i].real << setw(15) << radix2.x[i].imag << endl;
	}
	cout << " output radix-4: " << endl;
	for (int i = 0; i < POINTS; i++) {
		cout << i << fixed << setw(15) << radix4.x[i].real << setw(15) << radix4.x[i].imag << endl;
	}
	cout << " output radix-8: " << endl;
	for (int i = 0; i < POINTS; i++) {
		cout << i << fixed << setw(15) << radix8.x[i].real << setw(15) << radix8.x[i].imag << endl;
	}
	//cout << " output " << endl;
	//for (int i = 0; i < POINTS; i++) {
	//	cout << i << setw(15) << dft.out[i].real << " , " << setw(15) << radix4.x[i].real << " , " << setw(15) << dft.out[i].imag << setw(15) << radix4.x[i].imag << " , " << endl;
	//}

	//cout << " output " << endl;
	//for (int i = 0; i < POINTS; i++) {
	//	cout << i << setw(15) << dft.out[i].real << " , " << setw(15) << radix8.x[i].real << " , " << setw(15) << dft.out[i].imag << setw(15) << radix8.x[i].imag << " , " << endl;
	//}
	//for (int i = 0; i < POINTS; i++) {
	//	dft.x[i] = Complex<TYPE>(i, 0);
	//}
	//dft.bitReverse();
	//for (int i = 0; i < POINTS; i++) {
	//	cout << i << setw(15) << dft.x[i].real << endl;
	//}

    return 0;
}
