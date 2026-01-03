#pragma once

#include "../../basics.hpp"

namespace numtheo_n {
	template<i64 P, bool _64> i32 legendre(ModIntPr<P, _64> x) {
		if (x == 0) {
			return 0;
		}
		if (qpow(x, (ModIntPr<P, _64>::mod() - 1) >> 1) == 1) {
			return 1;
		}
		return -1;
	}
	template<i64 P, bool _64> std::optional<ModIntPr<P, _64>> sqrt(ModIntPr<P, _64> x) {
		if (ModIntPr<P, _64>::mod() == 2) {
			return x;
		}
		if (x == 0) {
			return 0;
		}
		if (legendre(x) == -1) {
			return std::nullopt;
		}
		std::mt19937_64 rndu(std::random_device{}());
		ModIntPr<P, _64> u, t;
		do {
			u = rndu();
			t = u * u - x;
		} while (legendre(t) != -1);
		struct sqrtt {
			ModIntPr<P, _64> re, im, _t;
			sqrtt &operator*=(sqrtt _x) {
				ModIntPr<P, _64> _r = re * _x.re + im * _x.im * _t;
				ModIntPr<P, _64> _i = re * _x.im + im * _x.re;
				re = _r, im = _i;
				return *this;
			}
		};
		ModIntPr<P, _64> ret = qpow(
			sqrtt{u, ModIntPr<P, _64>(1, false), t},
			(ModIntPr<P, _64>::mod() + 1) >> 1,
			sqrtt{ModIntPr<P, _64>(1, false), ModIntPr<P, _64>(0, false), t}
		).re;
		if (ret.val > (-ret).val) {
			return -ret;
		} else {
			return ret;
		}
	};
}