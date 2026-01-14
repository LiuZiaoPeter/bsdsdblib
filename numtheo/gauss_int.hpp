#pragma once

#include <algorithm>
#include <complex>
#include <iostream>

#include "../basics.hpp"

namespace numtheo {
	template<i128::signed_integral T> struct GaussInt {
		T re, im;
		GaussInt(T r = 0, T i = 0) : re(r), im(i) {}
		template<class U> operator std::complex<U>() const requires std::convertible_to<T, U> {
			return std::complex<U>(static_cast<U>(re), static_cast<U>(im));
		}
		void outp_onlyspace(std::ostream& out = std::cout) {
			out << re << ' ' << im;
		}
		void outp_bracecomma(std::ostream& out = std::cout) {
			out << '(' << re << ',' << im << ')';
		}
	};
	template<i128::signed_integral T> GaussInt<T> conj(const GaussInt<T> x) {
		return GaussInt<T>(x.re, -x.im);
	}
	template<i128::signed_integral T> GaussInt<T> operator+(const GaussInt<T> x, const GaussInt<T> y) {
		return GaussInt<T>(x.re + y.re, x.im + y.im);
	}
	template<i128::signed_integral T> GaussInt<T> operator-(const GaussInt<T> x) {
		return GaussInt<T>(-x.re, -x.im);
	}
	template<i128::signed_integral T> GaussInt<T> operator-(const GaussInt<T> x, const GaussInt<T> y) {
		return GaussInt<T>(x.re - y.re, x.im - y.im);
	}
	template<i128::signed_integral T> GaussInt<T> operator*(const GaussInt<T> x, const GaussInt<T> y) {
		return GaussInt<T>(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re);
	}
	template<i128::liftable_signed T> i128::up_t<T> norm(const GaussInt<T> x) {
		return static_cast<i128::up_t<T>>(x.re) * x.re + static_cast<i128::up_t<T>>(x.im) * x.im;
	}
	template<i128::signed_integral T> GaussInt<T> operator/(const GaussInt<T> x, const GaussInt<T> y) {
		std::complex<double> retc = static_cast<std::complex<double>>(x) / static_cast<std::complex<double>>(y);
		return GaussInt<T>(static_cast<T>(std::round(retc.real())), static_cast<T>(std::round(retc.imag())));
	}
	template<i128::signed_integral T> GaussInt<T> operator%(const GaussInt<T> x, const GaussInt<T> y) {
		return x - x / y * y;
	}
	template<i128::signed_integral T> GaussInt<T> &operator+=(GaussInt<T> &x, const GaussInt<T> y) {
		return x = x + y;
	}
	template<i128::signed_integral T> GaussInt<T> &operator-=(GaussInt<T> &x, const GaussInt<T> y) {
		return x = x - y;
	}
	template<i128::signed_integral T> GaussInt<T> &operator*=(GaussInt<T> &x, const GaussInt<T> y) {
		return x = x * y;
	}
	template<i128::signed_integral T> GaussInt<T> &operator/=(GaussInt<T> &x, const GaussInt<T> y) {
		return x = x / y;
	}
	template<i128::signed_integral T> GaussInt<T> &operator%=(GaussInt<T> &x, const GaussInt<T> y) {
		return x = x % y;
	}
	template<i128::signed_integral T> GaussInt<T> gauss_gcd(GaussInt<T> x, GaussInt<T> y) {
		if (norm(x) < norm(y)) {
			std::swap(x, y);
		}
		while (norm(y)) {
			x %= y;
			std::swap(x, y);
		}
		return x;
	}
}
