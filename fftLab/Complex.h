//

#pragma once

template <typename T>
class Complex {
public:
	Complex() : real(0), imag(0) {}
	template <typename U>
	Complex(const U& r) : real(static_cast<T>(r)), imag(static_cast<T>(0)) {}
	template <typename U, typename V>
	Complex(const U& r, const V&i) : real(static_cast<T>(r)), imag(static_cast<T>(i)) {}
	~Complex() {}

	__forceinline Complex<T> operator+ (const Complex<T>& b) const {
		return Complex(real + b.real, imag + b.imag);
	}
	__forceinline Complex<T> operator- (const Complex<T>& b) const {
		return Complex(real - b.real, imag - b.imag);
	}
	__forceinline Complex<T> operator* (const Complex<T>& b) const {
		return Complex<T>((real*b.real - imag*b.imag), (real*b.imag + imag*b.real));
	}
	__forceinline void operator+= (const Complex<T>& b) {
		real += b.real;	imag += b.imag;
	}
	__forceinline void operator-= (const Complex<T>& b) {
		real -= b.real;	imag -= b.imag;
	}

public:
	T real;
	T imag;
};

