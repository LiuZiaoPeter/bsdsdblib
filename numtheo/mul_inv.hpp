#pragma once

#include <cmath>
#include <type_traits>
#include <vector>

#include "../basics.hpp"
#include "../int128.hpp"
#include "modint.hpp"

namespace mul_inv_hpp {
	template<i64 P, bool _64> struct aux {
		using val_t = std::conditional_t<_64, u64, u32>;
		inline static std::vector<numtheo_n::ModIntPr<P, _64>> inv_v;
		inline static val_t cbrtP_log2, cbrtP, cbrtP2;
		inline static bool O1inv_mode = false;
		inline static std::vector<std::pair<val_t, val_t>> farey_v, farey_prec, farey_succ;
	};
}

namespace numtheo_n {
	template<i128::signed_integral T, i64 P, bool _64> ModIntPr<P, _64> qpow_signed(ModIntPr<P, _64> x, T y) {
		using MIP = ModIntPr<P, _64>;
		y %= MIP::mod() - 1;
		if (y < 0) {
			y += MIP::mod() - 1;
		}
		return qpow(x, y);
	}
	template<i64 P, bool _64, i128::unsigned_integral T> void lin_inv_preproc(T N) {
		using MIP = ModIntPr<P, _64>;
		using aux = mul_inv_hpp::aux<P, _64>;
		aux::inv_v.resize(N + 1);
		aux::inv_v[0] = 0, aux::inv_v[1] = 1;
		for (u32 i = 2; i <= N; ++i) {
			aux::inv_v[i] = (-MIP(MIP::mod() / i, false)) * aux::inv_v[MIP::mod() % i];
		}
	}
	template<i64 P, bool _64> void O1inv_preproc() {
		using MIP = ModIntPr<P, _64>;
		using aux = mul_inv_hpp::aux<P, _64>;
		using val_t = typename MIP::val_t;
		aux::O1inv_mode = true;
		aux::cbrtP_log2 = std::bit_width(static_cast<val_t>(std::cbrt(MIP::mod()) + 1));
		aux::cbrtP = static_cast<val_t>(1) << aux::cbrtP_log2;
		aux::cbrtP2 = static_cast<val_t>(1) << (aux::cbrtP_log2 << 1);
		if (aux::inv_v.size() <= aux::cbrtP2) {
			lin_inv_preproc<P, _64>(aux::cbrtP2);
		}
		aux::farey_v.assign(aux::cbrtP2 + 1, std::make_pair(0, 0));
		for (val_t p = 0; p <= aux::cbrtP; ++p) {
			for (val_t q = (p == 1 ? 1 : p + 1); q <= aux::cbrtP; ++q) {
				val_t cur = (p << (aux::cbrtP_log2 << 1)) / q;
				if (aux::farey_v[cur].second == 0) {
					aux::farey_v[cur] = std::make_pair(p, q);
				}
			}
		}
		aux::farey_prec = aux::farey_v;
		for (auto it = aux::farey_prec.begin(); it != aux::farey_prec.end(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
		aux::farey_succ = aux::farey_v;
		for (auto it = aux::farey_succ.rbegin(); it != aux::farey_succ.rend(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
	}
	template<i64 P, bool _64> ModIntPr<P, _64> inv(ModIntPr<P, _64> x) {
		using MIP = ModIntPr<P, _64>;
		using aux = mul_inv_hpp::aux<P, _64>;
		using val_t = typename MIP::val_t;
		using sval_t = i128::make_signed_t<val_t>;
		using mul_t = typename MIP::mul_t;
		using smul_t = i128::make_signed_t<mul_t>;
		if (x.value() == 0) {
			throw std::invalid_argument(
				"numtheo_n::template<i64 P, bool _64> ModIntPr<P, _64> inv(void) : getting inv(0)"
			);
		}
		if (x.value() < aux::inv_v.size()) {
			return aux::inv_v[x.value()];
		}
		if (aux::O1inv_mode == true) {
			val_t cur = static_cast<val_t>((static_cast<mul_t>(x.value()) << (aux::cbrtP_log2 << 1)) / MIP::mod());
			val_t q;
			sval_t t = static_cast<sval_t>(static_cast<smul_t>(x.value()) * aux::farey_prec[cur].second
					- static_cast<smul_t>(MIP::mod()) * aux::farey_prec[cur].first);
			if (static_cast<val_t>(std::abs(t)) <= aux::cbrtP2) {
				q = aux::farey_prec[cur].second;
			} else {
				q = aux::farey_succ[cur].second;
				t = static_cast<sval_t>(static_cast<smul_t>(x.value()) * q
					- static_cast<smul_t>(MIP::mod()) * aux::farey_succ[cur].first);
			}
			MIP tinv = ((t >= 0) ? aux::inv_v[t] : -aux::inv_v[-t]);
			return MIP(q, false) * tinv;
		}
		return qpow(x, MIP::mod() - 2);
	}
	template<i64 P, bool _64> ModIntPr<P, _64> operator/(const ModIntPr<P, _64> x, const ModIntPr<P, _64> y) {
		return x * inv(y);
	}
	template<i64 P, bool _64> ModIntPr<P, _64> &operator/=(ModIntPr<P, _64> &x, const ModIntPr<P, _64> y) {
		return x *= inv(y);
	}
}