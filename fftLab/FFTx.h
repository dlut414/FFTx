//

#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "../common/Complex.h"

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

	template <typename R>
	__forceinline C EI(const R& k) {
		return C(cos(k), sin(k));
	}
	__forceinline C W(const unsigned& k) const {
		const T PI = T(3.14159265358979323846);
		const T theta = static_cast<T>(-2 * PI / k);
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
		const T PI = T(3.14159265358979323846);
		for (unsigned i = 0; i < N; i++) {
			out[i] = C();
			for (unsigned j = 0; j < N; j++) {
				out[i] += x[j] * EI(-2*PI*i*j / N);
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

	void Run_DIT() {
		BitReverse(N, x);
		Perform_DIT(N, x);
	}
	void Run_DIT_aligned() {
		BitReverse(N, x);
		Perform_DIT_aligned(N, x);
	}
	void Run_DIF() {
		Perform_DIF(N, x);
		BitReverse(N, x);
	}
	void Run_DIF_aligned() {
		Perform_DIF_aligned(N, x);
		BitReverse(N, x);
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform_DIT(const unsigned& nn, C* const data) {
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
	void Perform_DIT_aligned(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const unsigned half = nn >> 1;
		std::vector<C> shift(nn);
		for (unsigned step = 1, groupSize = nn; step < (nn); step <<= 1, groupSize >>= 1) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step))* factorId;
				const C factor = EI(theta);
				const unsigned offset = factorId* groupSize;
				for (unsigned i = 0; i < groupSize; i+=2) {
					const unsigned a = offset + i;
					const unsigned b = a + 1;
					const unsigned aa = a >> 1;
					const C wb = factor*data[b];
					shift[aa] = data[a] + wb;
					shift[aa+half] = data[a] - wb;
				}
			}
			for (unsigned i = 0; i < nn; i++) {
				data[i] = shift[i];
			}
		}
	}
	void Perform_DIF(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		for (unsigned step = nn >> 1; step >= 1; step >>= 1) {
			const unsigned jump = step << 1;
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step))* factorId;
				const C factor = EI(theta);
				for (unsigned a = factorId; a<(nn); a += jump) {
					const unsigned b = a + step;
					const C res0 = data[a] + data[b];
					const C res1 = (data[a] - data[b])* factor;
					data[a] = res0;
					data[b] = res1;
				}
			}
		}
	}
	void Perform_DIF_aligned(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const unsigned half = nn >> 1;
		std::vector<C> shift(nn);
		for (unsigned step = nn >> 1, groupSize = 2; step >= 1; step >>= 1, groupSize <<= 1) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step))* factorId;
				const C factor = EI(theta);
				const unsigned offset = factorId* groupSize;
				for (unsigned i = 0; i < groupSize; i += 2) {
					const unsigned a = offset + i;
					const unsigned b = a + 1;
					const unsigned aa = a >> 1;
					shift[a] = (data[aa] + data[aa + half]);
					shift[b] = (data[aa] - data[aa + half])*factor;
				}
			}
			for (unsigned i = 0; i < nn; i++) {
				data[i] = shift[i];
			}
		}
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

	void Run_DIT() {
		BitReverse(N, x);
		Perform_DIT(N, x);
	}
	void Run_DIT_aligned() {
		BitReverse(N, x);
		Perform_DIT_aligned(N, x);
	}
	void Run_DIF() {
		Perform_DIF(N, x);
		BitReverse(N, x);
	}
	void Run_DIF_aligned() {
		Perform_DIF_aligned(N, x);
		BitReverse(N, x);
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform_DIT(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C jj = C(0, 1);
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
					C res0 = data[a] +		factor1*data[b] +	factor2*data[c] +		factor3*data[d];
					C res1 = data[a] - jj*	factor1*data[b] -	factor2*data[c] +	jj*	factor3*data[d];
					C res2 = data[a] -		factor1*data[b] +	factor2*data[c] -		factor3*data[d];
					C res3 = data[a] + jj*	factor1*data[b] -	factor2*data[c] -	jj*	factor3*data[d];
					data[a] = res0;
					data[b] = res1;
					data[c] = res2;
					data[d] = res3;
				}
			}
		}
	}
	void Perform_DIT_aligned(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C jj = C(0, 1);
		const unsigned quarter = nn >> 2;
		std::vector<C> shift(nn);
		for (unsigned step = 1, groupSize = nn; step < (nn); step <<= 2, groupSize >>= 2) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 1))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta * 2);
				const C factor3 = EI(theta * 3);
				const unsigned offset = factorId* groupSize;
				for (unsigned i = 0; i < groupSize; i += 4) {
					const unsigned a = offset + i;
					const unsigned b = a + 1;
					const unsigned c = a + 2;
					const unsigned d = a + 3;
					const unsigned aa = a >> 2;
					const unsigned bb = aa + quarter;
					const unsigned cc = aa + quarter * 2;
					const unsigned dd = aa + quarter * 3;
					shift[aa] = data[a] +			factor1*data[b] +			factor2*data[c] +			factor3*data[d];
					shift[bb] = data[a] - jj*		factor1*data[b] -			factor2*data[c] + jj*		factor3*data[d];
					shift[cc] = data[a] -			factor1*data[b] +			factor2*data[c] -			factor3*data[d];
					shift[dd] = data[a] + jj*		factor1*data[b] -			factor2*data[c] - jj*		factor3*data[d];
				}
			}
			for (unsigned i = 0; i < nn; i++) {
				data[i] = shift[i];
			}
		}
	}
	void Perform_DIF(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C jj = C(0, 1);
		for (unsigned step = nn >> 2; step >= 1; step >>= 2) {
			const unsigned jump = step << 2;
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 1))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta * 2);
				const C factor3 = EI(theta * 3);
				for (unsigned a = factorId; a<(nn); a += jump) {
					const unsigned b = a + step;
					const unsigned c = a + step * 2;
					const unsigned d = a + step * 3;
					C res0 = (data[a] +		data[b] +	data[c] +		data[d]);
					C res1 = (data[a] -	jj*	data[b] -	data[c] +	jj*	data[d])*factor1;
					C res2 = (data[a] -		data[b] +	data[c] -		data[d])*factor2;
					C res3 = (data[a] + jj*	data[b] -	data[c] -	jj*	data[d])*factor3;
					data[a] = res0;
					data[b] = res1;
					data[c] = res2;
					data[d] = res3;
				}
			}
		}
	}
	void Perform_DIF_aligned(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C jj = C(0, 1);
		const unsigned quarter = nn >> 2;
		std::vector<C> shift(nn);
		for (unsigned step = nn >> 2, groupSize = 4; step >= 1; step >>= 2, groupSize <<= 2) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 1))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta * 2);
				const C factor3 = EI(theta * 3);
				const unsigned offset = factorId* groupSize;
				for (unsigned i = 0; i < groupSize; i += 4) {
					const unsigned a = offset + i;
					const unsigned b = a + 1;
					const unsigned c = a + 2;
					const unsigned d = a + 3;
					const unsigned aa = a >> 2;
					const unsigned bb = aa + quarter;
					const unsigned cc = aa + quarter * 2;
					const unsigned dd = aa + quarter * 3;
					shift[a] = (data[aa] +			data[bb] +			data[cc] +			data[dd]);
					shift[b] = (data[aa] - jj*		data[bb] -			data[cc] + jj*		data[dd])*factor1;
					shift[c] = (data[aa] -			data[bb] +			data[cc] -			data[dd])*factor2;
					shift[d] = (data[aa] + jj*		data[bb] -			data[cc] - jj*		data[dd])*factor3;
				}
			}
			for (unsigned i = 0; i < nn; i++) {
				data[i] = shift[i];
			}
		}
	}
	__forceinline void initialize() {
		int flag = 0;
		std::vector<unsigned> tempt(N);
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
			if (mirror > i) {
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

	void Run_DIT() {
		BitReverse(N, x);
		Perform_DIT(N, x);
	}
	void Run_DIT_aligned() {
		BitReverse(N, x);
		Perform_DIT_aligned(N, x);
	}
	void Run_DIF() {
		Perform_DIF(N, x);
		BitReverse(N, x);
	}
	void Run_DIF_aligned() {
		Perform_DIF_aligned(N, x);
		BitReverse(N, x);
	}

protected:
	void Perform(const unsigned& nn, const C* const input, C* const output) {

	}
	void Perform_DIT(const unsigned& nn, C* const data) {
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
	void Perform_DIT_aligned(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C W8 = W(8);
		const C jW8 = C(0, 1)*W8;
		const C jj = C(0, 1);
		const unsigned fraction = nn >> 3;
		std::vector<C> shift(nn);
		for (unsigned step = 1, groupSize = nn; step < (nn); step <<= 3, groupSize >>= 3) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 2))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta * 2);
				const C factor3 = EI(theta * 3);
				const C factor4 = EI(theta * 4);
				const C factor5 = EI(theta * 5);
				const C factor6 = EI(theta * 6);
				const C factor7 = EI(theta * 7);
				const unsigned offset = factorId* groupSize;
				for (unsigned i = 0; i < groupSize; i += 8) {
					const unsigned a0 = offset + i;
					const unsigned a1 = a0 + 1;
					const unsigned a2 = a0 + 2;
					const unsigned a3 = a0 + 3;
					const unsigned a4 = a0 + 4;
					const unsigned a5 = a0 + 5;
					const unsigned a6 = a0 + 6;
					const unsigned a7 = a0 + 7;
					const unsigned aa0 = a0 >> 3;
					const unsigned aa1 = aa0 + fraction;
					const unsigned aa2 = aa0 + fraction * 2;
					const unsigned aa3 = aa0 + fraction * 3;
					const unsigned aa4 = aa0 + fraction * 4;
					const unsigned aa5 = aa0 + fraction * 5;
					const unsigned aa6 = aa0 + fraction * 6;
					const unsigned aa7 = aa0 + fraction * 7;
					shift[aa0] = data[a0] +		factor1*data[a1] +		factor2*data[a2] +		factor3*data[a3] + factor4*data[a4] +		factor5*data[a5] +		factor6*data[a6] +		factor7*data[a7];
					shift[aa1] = data[a0] + W8*	factor1*data[a1] - jj*	factor2*data[a2] - jW8*	factor3*data[a3] - factor4*data[a4] - W8*	factor5*data[a5] + jj*	factor6*data[a6] + jW8*	factor7*data[a7];
					shift[aa2] = data[a0] - jj*	factor1*data[a1] -		factor2*data[a2] + jj*	factor3*data[a3] + factor4*data[a4] - jj*	factor5*data[a5] -		factor6*data[a6] + jj*	factor7*data[a7];
					shift[aa3] = data[a0] - jW8*factor1*data[a1] + jj*	factor2*data[a2] + W8*	factor3*data[a3] - factor4*data[a4] + jW8*	factor5*data[a5] - jj*	factor6*data[a6] - W8*	factor7*data[a7];
					shift[aa4] = data[a0] -		factor1*data[a1] +		factor2*data[a2] -		factor3*data[a3] + factor4*data[a4] -		factor5*data[a5] +		factor6*data[a6] -		factor7*data[a7];
					shift[aa5] = data[a0] - W8*	factor1*data[a1] - jj*	factor2*data[a2] + jW8*	factor3*data[a3] - factor4*data[a4] + W8*	factor5*data[a5] + jj*	factor6*data[a6] - jW8*	factor7*data[a7];
					shift[aa6] = data[a0] + jj*	factor1*data[a1] -		factor2*data[a2] - jj*	factor3*data[a3] + factor4*data[a4] + jj*	factor5*data[a5] -		factor6*data[a6] - jj*	factor7*data[a7];
					shift[aa7] = data[a0] + jW8*factor1*data[a1] + jj*	factor2*data[a2] - W8*	factor3*data[a3] - factor4*data[a4] - jW8*	factor5*data[a5] - jj*	factor6*data[a6] + W8*	factor7*data[a7];
				}
			}
			for (unsigned i = 0; i < nn; i++) {
				data[i] = shift[i];
			}
		}
	}
	void Perform_DIF(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C W8 = W(8);
		const C jW8 = C(0, 1)*W8;
		const C jj = C(0, 1);
		for (unsigned step = nn >> 3; step >= 1; step >>= 3) {
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
					C res0 = (data[a0] +		data[a1] +		data[a2] +		data[a3] + data[a4] +		data[a5] +		data[a6] +		data[a7]);
					C res1 = (data[a0] + W8*	data[a1] - jj*	data[a2] - jW8*	data[a3] - data[a4] - W8*	data[a5] + jj*	data[a6] + jW8*	data[a7])*factor1;
					C res2 = (data[a0] - jj*	data[a1] -		data[a2] + jj*	data[a3] + data[a4] - jj*	data[a5] -		data[a6] + jj*	data[a7])*factor2;
					C res3 = (data[a0] - jW8*	data[a1] + jj*	data[a2] + W8*	data[a3] - data[a4] + jW8*	data[a5] - jj*	data[a6] - W8*	data[a7])*factor3;
					C res4 = (data[a0] -		data[a1] +		data[a2] -		data[a3] + data[a4] -		data[a5] +		data[a6] -		data[a7])*factor4;
					C res5 = (data[a0] - W8*	data[a1] - jj*	data[a2] + jW8*	data[a3] - data[a4] + W8*	data[a5] + jj*	data[a6] - jW8*	data[a7])*factor5;
					C res6 = (data[a0] + jj*	data[a1] -		data[a2] - jj*	data[a3] + data[a4] + jj*	data[a5] -		data[a6] - jj*	data[a7])*factor6;
					C res7 = (data[a0] + jW8*	data[a1] + jj*	data[a2] - W8*	data[a3] - data[a4] - jW8*	data[a5] - jj*	data[a6] + W8*	data[a7])*factor7;
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
	void Perform_DIF_aligned(const unsigned& nn, C* const data) {
		const T PI = T(3.14159265358979323846);
		const C W8 = W(8);
		const C jW8 = C(0, 1)*W8;
		const C jj = C(0, 1);
		const unsigned fraction = nn >> 3;
		std::vector<C> shift(nn);
		for (unsigned step = nn >> 3, groupSize = 8; step >= 1; step >>= 3, groupSize <<= 3) {
			for (unsigned factorId = 0; factorId < step; factorId++) {
				const T theta = -(PI / T(step << 2))* factorId;
				const C factor1 = EI(theta);
				const C factor2 = EI(theta * 2);
				const C factor3 = EI(theta * 3);
				const C factor4 = EI(theta * 4);
				const C factor5 = EI(theta * 5);
				const C factor6 = EI(theta * 6);
				const C factor7 = EI(theta * 7);
				const unsigned offset = factorId* groupSize;
				for (unsigned i = 0; i < groupSize; i += 8) {
					const unsigned a0 = offset + i;
					const unsigned a1 = a0 + 1;
					const unsigned a2 = a0 + 2;
					const unsigned a3 = a0 + 3;
					const unsigned a4 = a0 + 4;
					const unsigned a5 = a0 + 5;
					const unsigned a6 = a0 + 6;
					const unsigned a7 = a0 + 7;
					const unsigned aa0 = a0 >> 3;
					const unsigned aa1 = aa0 + fraction;
					const unsigned aa2 = aa0 + fraction * 2;
					const unsigned aa3 = aa0 + fraction * 3;
					const unsigned aa4 = aa0 + fraction * 4;
					const unsigned aa5 = aa0 + fraction * 5;
					const unsigned aa6 = aa0 + fraction * 6;
					const unsigned aa7 = aa0 + fraction * 7;
					shift[a0] = (data[aa0] +		data[aa1] +		data[aa2] +			data[aa3] + data[aa4] +			data[aa5] +		data[aa6] +			data[aa7]);
					shift[a1] = (data[aa0] + W8*	data[aa1] - jj*	data[aa2] - jW8*	data[aa3] - data[aa4] - W8*		data[aa5] + jj*	data[aa6] + jW8*	data[aa7])*factor1;
					shift[a2] = (data[aa0] - jj*	data[aa1] -		data[aa2] + jj*		data[aa3] + data[aa4] - jj*		data[aa5] -		data[aa6] + jj*		data[aa7])*factor2;
					shift[a3] = (data[aa0] - jW8*	data[aa1] + jj*	data[aa2] + W8*		data[aa3] - data[aa4] + jW8*	data[aa5] - jj*	data[aa6] - W8*		data[aa7])*factor3;
					shift[a4] = (data[aa0] -		data[aa1] +		data[aa2] -			data[aa3] + data[aa4] -			data[aa5] +		data[aa6] -			data[aa7])*factor4;
					shift[a5] = (data[aa0] - W8*	data[aa1] - jj*	data[aa2] + jW8*	data[aa3] - data[aa4] + W8*		data[aa5] + jj*	data[aa6] - jW8*	data[aa7])*factor5;
					shift[a6] = (data[aa0] + jj*	data[aa1] -		data[aa2] - jj*		data[aa3] + data[aa4] + jj*		data[aa5] -		data[aa6] - jj*		data[aa7])*factor6;
					shift[a7] = (data[aa0] + jW8*	data[aa1] + jj*	data[aa2] - W8*		data[aa3] - data[aa4] - jW8*	data[aa5] - jj*	data[aa6] + W8*		data[aa7])*factor7;
				}
			}
			for (unsigned i = 0; i < nn; i++) {
				data[i] = shift[i];
			}
		}
	}
	__forceinline void initialize() {
		int flag = 0;
		std::vector<unsigned> tempt(N);
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