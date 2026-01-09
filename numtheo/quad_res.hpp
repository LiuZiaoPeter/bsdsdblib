// tested by lg_P5491

#pragma once

#include <random>
#include <optional>

#include "../basics.hpp"
#include "modint.hpp"

namespace numtheo {
	template<i64 P, bool _64> i32 legendre(ModIntPr<P, _64> x) {
		using MIP = ModIntPr<P, _64>;
		if (x == 0) {
			return 0;
		}
		if (qpow(x, (MIP::mod() - 1) >> 1, MIP(1, false)) == 1) {
			return 1;
		}
		return -1;
	}
	template<i64 P, bool _64> std::optional<ModIntPr<P, _64>> sqrt(ModIntPr<P, _64> x) {
		using MIP = ModIntPr<P, _64>;
		if (MIP::mod() == 2) {
			return x;
		}
		if (x == 0) {
			return 0;
		}
		if (legendre(x) == -1) {
			return std::nullopt;
		}
		std::mt19937_64 rndu(std::random_device{}());
		MIP u, t;
		do {
			u = rndu();
			t = u * u - x;
		} while (legendre(t) != -1);
		struct sqrtt {
			MIP re, im, _t;
			sqrtt &operator*=(sqrtt _x) {
				MIP _r = re * _x.re + im * _x.im * _t;
				MIP _i = re * _x.im + im * _x.re;
				re = _r, im = _i;
				return *this;
			}
		};
		MIP ret = qpow(sqrtt{u, MIP(1, false), t}, (MIP::mod() + 1) >> 1, sqrtt{MIP(1, false), MIP(0, false), t}).re;
		if (ret.value() > (MIP::mod() >> 1)) {
			return -ret;
		} else {
			return ret;
		}
	};
}