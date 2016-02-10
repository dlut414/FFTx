//

#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "Complex.h"

template <unsigned BIT>	struct POWER2		{ enum { value = 2 * POWER2<BIT-1>::value, }; };
template <>			struct POWER2<1>	{ enum { value = 2, }; };
template <unsigned BIT>	struct POWER4		{ enum { value = 4 * POWER4<BIT-1>::value, }; };
template <>			struct POWER4<1>	{ enum { value = 4, }; };
template <unsigned BIT>	struct POWER8		{ enum { value = 8 * POWER8<BIT-1>::value, }; };
template <>			struct POWER8<1>	{ enum { value = 8, }; };

template <typename T, unsigned BIT, unsigned N>
class FFTx {
public:
	struct Points { enum { value = N, }; };
	typedef Complex<T> C;
	FFTx() { initialize(); }
	~FFTx() {}

	virtual void Run() = 0;

	template <typename R>
	__forceinline C EI(const R& k) {
		return C(cos(k), sin(k));
	}
	__forceinline C W(const unsigned& k) const {
		const T theta = static_cast<T>(-2 * M_PI* k / N);
		return C(cos(theta), sin(theta));
	}
	__forceinline void BitReverse(const unsigned& nn, C* const data) {
		for (unsigned i = 0; i < nn; i++) {
			unsigned mirror = BitMirror[i];
			if (mirror > unsigned(i)) {
				C temp = data[i];
				data[i] = data[mirror];
				data[mirror] = temp;
			}
		}
	}
protected:
	void initialize() {
		unsigned mirror = 0;
		for (unsigned i = 0; i < N; i++) {
			BitMirror[i] = mirror;
			unsigned mask = N;
			while (mirror & (mask >>= 1)) mirror &= ~mask;
			mirror |= mask;
		}
	}

public:
	C x[N];
	unsigned BitMirror[N];
};

template <typename T, unsigned BIT, unsigned N = POWER2<BIT>::value>
class DFT : public FFTx<T,BIT,N> {
public:
	DFT() {}
	~DFT() {}

	void Run() {
		for (unsigned i = 0; i < N; i++) {
			out[i] = C();
			for (unsigned j = 0; j < N; j++) {
				out[i] += x[j] * W(i*j);
			}
		}
	}

public:
	C out[N];

};

template <typename T, unsigned BIT, unsigned N = POWER2<BIT>::value>
class Radix_2 : public FFTx<T,BIT,N> {
public:
	Radix_2() {}
	~Radix_2() {}

	void Run() {
		BitReverse(N, x);
		Perform(N, x);
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		for (unsigned step = 1; step < (nn); step <<= 1) {
			const unsigned jump = step << 1;
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const C factor = EI(-(PI / T(step))* factorId);
				for (unsigned a = factorId; a<(nn); a += jump) {
					const unsigned b = a + step;
					const C wb = factor*data[b];
					data[b] = data[a] - wb;
					data[a] = data[a] + wb;
				}
			}
		}
	}
	void Perform_aligned(const unsigned& nn, const C* const input, C* const output) {
		const T PI = T(3.14159265358979323846);
		const unsigned jump = 1 << 1;
		const unsigned half = nn >> 1;
		for (unsigned step = 1; step < (nn); step <<= 1) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const C factor = EI(-(PI / T(step))* factorId);
				for (unsigned a = 0; a<(nn); a += jump) {
					const unsigned b = a + 1;
					const unsigned f = a >> 1;
					const unsigned g = f + half;
					const C wb = factor* input[b];
					output[f] = input[a] - wb;
					output[g] = input[a] + wb;
				}
			}
		}
	}
};

template <typename T, unsigned BIT, unsigned N = POWER4<BIT>::value>
class Radix_4 : public FFTx<T,BIT,N> {
public:
	Radix_4() {}
	~Radix_4() {}

	void Run() {
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		for (unsigned step = 1; step < (nn); step <<= 2) {
			const unsigned jump = step << 2;
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const C factor1 = EI(-(PI / T(2*step))* factorId);
				const C factor2 = factor1*C(2);
				const C factor3 = factor1*C(3);
				for (unsigned a = factorId; a<(nn); a += jump) {
					const unsigned b = a + step;
					const unsigned c = a + step*2;
					const unsigned d = a + step*3;
					C res1 = data[a] + factor1*data[b] + factor2*data[c] + factor3*data[d];
					C res2 = data[a] + C(0, -1)*factor1*data[b] + C(-1, 0)*factor2*data[c] + C(0, 1)*factor3*data[d];
					C res3 = data[a] + C(-1, 0)*factor1*data[b] + C(1, 0)*factor2*data[c] + C(-1, 0)*factor3*data[d];
					C res4 = data[a] + C(0, 1)*factor1*data[b] + C(-1, 0)*factor2*data[c] + C(0, -1)*factor3*data[d];
					data[a] = res1;
					data[b] = res2;
					data[c] = res3;
					data[d] = res4;
				}
			}
		}
	}
};

template <typename T, unsigned BIT, unsigned N = POWER8<BIT>::value>
class Radix_8 : public FFTx<T,BIT,N> {
public:
	Radix_8() {}
	~Radix_8() {}

	void Run() {}
};