#pragma once

#include <cmath>
#include <utility>
#include <vector>

#include "../basics.hpp"
#include "modint.hpp"

namespace mul_inv_hpp {
	template<i64 P, bool _64> struct aux {
		inline static std::vector<numtheo::ModIntPr<P, _64>> inv_v;
		inline static u32 cbrtP_log2, cbrtP, cbrtP2;
		inline static bool O1inv_mode = false;
		inline static std::vector<std::pair<u32, u32>> farey_v, farey_prec, farey_succ;
	};
}

namespace numtheo {
	template<i64 P, bool _64, i128::signed_integral T> ModIntPr<P, _64> qpow_signed(ModIntPr<P, _64> x, T y) {
		using MIP = ModIntPr<P, _64>;
		y %= MIP::mod() - 1;
		if (y < 0) {
			y += MIP::mod() - 1;
		}
		return qpow(x, y, MIP(1, false));
	}
	template<i64 P, bool _64> void lin_inv_preproc(u32 N) {
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
		aux::O1inv_mode = true;
		aux::cbrtP_log2 = std::bit_width(static_cast<u32>(std::cbrt(MIP::mod()) + 1));
		aux::cbrtP = 1 << aux::cbrtP_log2;
		aux::cbrtP2 = 1 << (aux::cbrtP_log2 << 1);
		if (aux::inv_v.size() <= aux::cbrtP2) {
			lin_inv_preproc<P, _64>(aux::cbrtP2);
		}
		aux::farey_v.assign(aux::cbrtP2 + 1, std::make_pair(0, 0));
		for (u32 p = 0; p <= aux::cbrtP; ++p) {
			for (u32 q = (p == 1 ? 1 : p + 1); q <= aux::cbrtP; ++q) {
				u32 cur = static_cast<u32>((static_cast<u64>(p) << (aux::cbrtP_log2 << 1)) / q);
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
		if (x.value() == 0) {
			throw std::invalid_argument(
				__func_str__ + " : getting inv(0)"
			);
		}
		if (x.value() < aux::inv_v.size()) {
			return aux::inv_v[x.value()];
		}
		if (aux::O1inv_mode == true) {
			u32 cur = static_cast<u32>((static_cast<u64>(x.value()) << (aux::cbrtP_log2 << 1)) / MIP::mod());
			u32 q;
			i32 t = static_cast<i32>(static_cast<i64>(x.value()) * aux::farey_prec[cur].second
					- static_cast<i64>(MIP::mod()) * aux::farey_prec[cur].first);
			if (static_cast<u32>(std::abs(t)) <= aux::cbrtP2) {
				q = aux::farey_prec[cur].second;
			} else {
				q = aux::farey_succ[cur].second;
				t = static_cast<i32>(static_cast<i64>(x.value()) * q
					- static_cast<i64>(MIP::mod()) * aux::farey_succ[cur].first);
			}
			MIP tinv = ((t >= 0) ? aux::inv_v[t] : -aux::inv_v[-t]);
			return MIP(q, false) * tinv;
		}
		return qpow(x, MIP::mod() - 2, MIP(1, false));
	}
	template<i64 P, bool _64> ModIntPr<P, _64> operator/(const ModIntPr<P, _64> x, const ModIntPr<P, _64> y) {
		return x * inv(y);
	}
	template<i64 P, bool _64> ModIntPr<P, _64> &operator/=(ModIntPr<P, _64> &x, const ModIntPr<P, _64> y) {
		return x *= inv(y);
	}
}
