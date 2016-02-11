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
	FFTx() {}
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

public:
	C x[N];
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
	Radix_2() { initialize(); }
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
				const T theta = -(PI / T(step))* factorId;
				const C factor = EI(theta);
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

	}
	__forceinline void initialize() {
		unsigned mirror = 0;
		for (unsigned i = 0; i < N; i++) {
			BitMirror[i] = mirror;
			unsigned mask = N;
			while (mirror & (mask >>= 1)) mirror &= ~mask;
			mirror |= mask;
		}
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
	unsigned BitMirror[N];

};

template <typename T, unsigned BIT, unsigned N = POWER2<BIT>::value>
class Radix_4 : public FFTx<T,BIT,N> {
public:
	Radix_4() { initialize(); }
	~Radix_4() {}

	void Run() {
		BitReverse(N, x);
		Perform(N, x);
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		for (unsigned step = 1; step < (nn); step <<= 2) {
			const unsigned jump = step << 2;
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 1))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta*2);
				const C factor3 = EI(theta*3);
				for (unsigned a = factorId; a<(nn); a += jump) {
					const unsigned b = a + step;
					const unsigned c = a + step*2;
					const unsigned d = a + step*3;
					C res1 = data[a] +			factor1*data[b] +			factor2*data[c] +			factor3*data[d];
					C res2 = data[a] + C(0, -1)*factor1*data[b] + C(-1, 0)*	factor2*data[c] + C(0, 1)*	factor3*data[d];
					C res3 = data[a] + C(-1, 0)*factor1*data[b] + C(1, 0)*	factor2*data[c] + C(-1, 0)*	factor3*data[d];
					C res4 = data[a] + C(0, 1)*	factor1*data[b] + C(-1, 0)*	factor2*data[c] + C(0, -1)*	factor3*data[d];
					data[a] = res1;
					data[b] = res2;
					data[c] = res3;
					data[d] = res4;
				}
			}
		}
	}
	void Perform_aligned(const unsigned& nn, C* const data) {

	}
	__forceinline void initialize() {
		int flag = 0;
		unsigned tempt[N];
		for (unsigned i = 0; i < N; i++) BitMirror[i] = i;
		for (unsigned step = 1, t = N; step < N; step <<= 2, t >>= 2) {
			unsigned jump = t >> 2;
			flag = !flag;
			for (unsigned i = 0; i < step; i++) {
				for (unsigned j = 0; j < t; j++) {
					unsigned group = j % 4;
					unsigned offset = j / 4;
					unsigned id = group * jump + offset;
					if (flag) tempt[i*t + id] = BitMirror[i*t + j];
					else BitMirror[i*t + id] = tempt[i*t + j];
				}
			}
		}
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
	unsigned BitMirror[N];
};

template <typename T, unsigned BIT, unsigned N = POWER2<BIT>::value>
class Radix_8 : public FFTx<T,BIT,N> {
public:
	Radix_8() { initialize(); }
	~Radix_8() {}

	void Run() {
		BitReverse(N, x);
		Perform(N, x);
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C W8 = W(8);
		const C jW8 = C(0, 1)*W8;
		const C jj = C(0, 1);
		for (unsigned step = 1; step < (nn); step <<= 3) {
			const unsigned jump = step << 3;
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 2))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta * 2);
				const C factor3 = EI(theta * 3);
				const C factor4 = EI(theta * 4);
				const C factor5 = EI(theta * 5);
				const C factor6 = EI(theta * 6);
				const C factor7 = EI(theta * 7);
				for (unsigned a0 = factorId; a0<(nn); a0 += jump) {
					const unsigned a1 = a0 + step;
					const unsigned a2 = a0 + step * 2;
					const unsigned a3 = a0 + step * 3;
					const unsigned a4 = a0 + step * 4;
					const unsigned a5 = a0 + step * 5;
					const unsigned a6 = a0 + step * 6;
					const unsigned a7 = a0 + step * 7;
					C res0 = data[a0] +		factor1*data[a1] +		factor2*data[a2] +		factor3*data[a3] + factor4*data[a4] +		factor5*data[a5] +		factor6*data[a6] +		factor7*data[a7];
					C res1 = data[a0] + W8*	factor1*data[a1] - jj*	factor2*data[a2] - jW8*	factor3*data[a3] - factor4*data[a4] - W8*	factor5*data[a5] + jj*	factor6*data[a6] + jW8*	factor7*data[a7];
					C res2 = data[a0] - jj*	factor1*data[a1] -		factor2*data[a2] + jj*	factor3*data[a3] + factor4*data[a4] - jj*	factor5*data[a5] -		factor6*data[a6] + jj*	factor7*data[a7];
					C res3 = data[a0] - jW8*factor1*data[a1] + jj*	factor2*data[a2] + W8*	factor3*data[a3] - factor4*data[a4] + jW8*	factor5*data[a5] - jj*	factor6*data[a6] - W8*	factor7*data[a7];
					C res4 = data[a0] -		factor1*data[a1] +		factor2*data[a2] -		factor3*data[a3] + factor4*data[a4] -		factor5*data[a5] +		factor6*data[a6] -		factor7*data[a7];
					C res5 = data[a0] - W8*	factor1*data[a1] - jj*	factor2*data[a2] + jW8*	factor3*data[a3] - factor4*data[a4] + W8*	factor5*data[a5] + jj*	factor6*data[a6] - jW8*	factor7*data[a7];
					C res6 = data[a0] + jj*	factor1*data[a1] -		factor2*data[a2] - jj*	factor3*data[a3] + factor4*data[a4] + jj*	factor5*data[a5] -		factor6*data[a6] - jj*	factor7*data[a7];
					C res7 = data[a0] + jW8*factor1*data[a1] + jj*	factor2*data[a2] - W8*	factor3*data[a3] - factor4*data[a4] - jW8*	factor5*data[a5] - jj*	factor6*data[a6] + W8*	factor7*data[a7];
					data[a0] = res0;
					data[a1] = res1;
					data[a2] = res2;
					data[a3] = res3;
					data[a4] = res4;
					data[a5] = res5;
					data[a6] = res6;
					data[a7] = res7;
				}
			}
		}
	}
	void Perform_aligned(const unsigned& nn, C* const data) {

	}
	__forceinline void initialize() {
		int flag = 0;
		unsigned tempt[N];
		for (unsigned i = 0; i < N; i++) BitMirror[i] = i;
		for (unsigned step = 1, t = N; step < N; step <<= 3, t >>= 3) {
			unsigned jump = t >> 3;
			flag = !flag;
			for (unsigned i = 0; i < step; i++) {
				for (unsigned j = 0; j < t; j++) {
					unsigned group = j % 8;
					unsigned offset = j / 8;
					unsigned id = group * jump + offset;
					if (flag) tempt[i*t + id] = BitMirror[i*t + j];
					else BitMirror[i*t + id] = tempt[i*t + j];
				}
			}
		}
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
	unsigned BitMirror[N];
};