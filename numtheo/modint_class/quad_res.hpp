#pragma once

#include "../../basics.hpp"

namespace numtheo_n {
	template<i32 P> i32 legendre(MIP<P> x) {
		if (x == 0) {
			return 0;
		}
		if (qpow(x, (MIP<P>::mod() - 1) >> 1) == 1) {
			return 1;
		}
		return -1;
	}
	template<i32 P> std::optional<MIP<P>> sqrt(MIP<P> x) {
		if (MIP<P>::mod() == 2) {
			return x;
		}
		if (x == 0) {
			return 0;
		}
		if (legendre(x) == -1) {
			return std::nullopt;
		}
		std::mt19937_64 rndu(std::random_device{}());
		MIP<P> u, t;
		do {
			u = rndu();
			t = u * u - x;
		} while (legendre(t) != -1);
		struct sqrtt {
			MIP<P> re, im, _t;
			sqrtt &operator*=(sqrtt _x) {
				MIP<P> _r = re * _x.re + im * _x.im * _t;
				MIP<P> _i = re * _x.im + im * _x.re;
				re = _r, im = _i;
				return *this;
			}
		};
		MIP<P> ret = qpow(
			sqrtt{u, MIP<P>(1, false), t},
			(MIP<P>::mod() + 1) >> 1,
			sqrtt{MIP<P>(1, false), MIP<P>(0, false), t}
		).re;
		if (ret.val > (-ret).val) {
			return -ret;
		} else {
			return ret;
		}
	};
}