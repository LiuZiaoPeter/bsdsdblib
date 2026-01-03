#pragma once

#include <vector>

#include "../../basics.hpp"
#include "../../int128.hpp"

namespace numtheo_n {
	template<i128::signed_integral T, i64 P, bool _64> ModIntPr<P, _64> qpow_signed(ModIntPr<P, _64> x, T y) {
		y %= ModIntPr<P, _64>::mod() - 1;
		if (y < 0) {
			y += ModIntPr<P, _64>::mod() - 1;
		}
		return qpow(x, y);
	}
	template<i64 P, bool _64> template<i128::unsigned_integral T> void ModIntPr<P, _64>::lin_inv_preproc(T N) {
		ModIntPr<P, _64>::inv_v.resize(N + 1);
		ModIntPr<P, _64>::inv_v[0] = 0, ModIntPr<P, _64>::inv_v[1] = 1;
		for (u32 i = 2; i <= N; ++i) {
			ModIntPr<P, _64>::inv_v[i] = (-ModIntPr<P, _64>(ModIntPr<P, _64>::mod() / i, false)) * ModIntPr<P, _64>::inv_v[ModIntPr<P, _64>::mod() % i];
		}
	}
	template<i64 P, bool _64> void ModIntPr<P, _64>::O1inv_preproc() {
		ModIntPr<P, _64>::O1inv_mode = true;
		ModIntPr<P, _64>::cbrtP_log2 = std::bit_width(
			static_cast<ModIntPr<P, _64>::val_t>(std::cbrt(ModIntPr<P, _64>::mod()) + 1)
		);
		ModIntPr<P, _64>::cbrtP = static_cast<ModIntPr<P, _64>::val_t>(1) << ModIntPr<P, _64>::cbrtP_log2;
		ModIntPr<P, _64>::cbrtP2 = static_cast<ModIntPr<P, _64>::val_t>(1) << (ModIntPr<P, _64>::cbrtP_log2 << 1);
		if (ModIntPr<P, _64>::inv_v.size() <= ModIntPr<P, _64>::cbrtP2) {
			lin_inv_preproc(ModIntPr<P, _64>::cbrtP2);
		}
		ModIntPr<P, _64>::farey_v.assign(ModIntPr<P, _64>::cbrtP2 + 1, std::make_pair(0, 0));
		for (ModIntPr<P, _64>::val_t p = 0; p <= ModIntPr<P, _64>::cbrtP; ++p) {
			for (ModIntPr<P, _64>::val_t q = (p == 1 ? 1 : p + 1); q <= ModIntPr<P, _64>::cbrtP; ++q) {
				ModIntPr<P, _64>::val_t cur = (p << (ModIntPr<P, _64>::cbrtP_log2 << 1)) / q;
				if (ModIntPr<P, _64>::farey_v[cur].second == 0) {
					ModIntPr<P, _64>::farey_v[cur] = std::make_pair(p, q);
				}
			}
		}
		ModIntPr<P, _64>::farey_prec = ModIntPr<P, _64>::farey_v;
		for (auto it = ModIntPr<P, _64>::farey_prec.begin(); it != ModIntPr<P, _64>::farey_prec.end(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
		ModIntPr<P, _64>::farey_succ = ModIntPr<P, _64>::farey_v;
		for (auto it = ModIntPr<P, _64>::farey_succ.rbegin(); it != ModIntPr<P, _64>::farey_succ.rend(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
	}
	template<i64 P, bool _64> ModIntPr<P, _64> inv(ModIntPr<P, _64> x) {
		if (x.val == 0) {
			throw std::invalid_argument(
				"numtheo_n::template<i32 P> ModIntPr<P, _64> inv(void) : getting inv(0)"
			);
		}
		if (x.val < ModIntPr<P, _64>::inv_v.size()) {
			return ModIntPr<P, _64>::inv_v[x.val];
		}
		if (ModIntPr<P, _64>::O1inv_mode == true) {
			typename ModIntPr<P, _64>::val_t cur = static_cast<ModIntPr<P, _64>::val_t>((static_cast<ModIntPr<P, _64>::mul_t>(x.val) << (ModIntPr<P, _64>::cbrtP_log2 << 1)) / ModIntPr<P, _64>::mod());
			typename ModIntPr<P, _64>::val_t q;
			i128::make_signed_t<typename ModIntPr<P, _64>::val_t> t = static_cast<i128::make_signed_t<typename ModIntPr<P, _64>::val_t>>(static_cast<i128::make_signed_t<typename ModIntPr<P, _64>::mul_t>>(x.val) * ModIntPr<P, _64>::farey_prec[cur].second - static_cast<i128::make_signed_t<typename ModIntPr<P, _64>::mul_t>>(ModIntPr<P, _64>::mod()) * ModIntPr<P, _64>::farey_prec[cur].first);
			if (static_cast<typename ModIntPr<P, _64>::val_t>(std::abs(t)) <= ModIntPr<P, _64>::cbrtP2) {
				q = ModIntPr<P, _64>::farey_prec[cur].second;
			} else {
				q = ModIntPr<P, _64>::farey_succ[cur].second;
				t = static_cast<i128::make_signed_t<typename ModIntPr<P, _64>::val_t>>(static_cast<i128::make_signed_t<typename ModIntPr<P, _64>::mul_t>>(x.val) * q - static_cast<i128::make_signed_t<typename ModIntPr<P, _64>::mul_t>>(ModIntPr<P, _64>::mod()) * ModIntPr<P, _64>::farey_succ[cur].first);
			}
			ModIntPr<P, _64> tinv = ((t >= 0) ? ModIntPr<P, _64>::inv_v[t] : -ModIntPr<P, _64>::inv_v[-t]);
			return ModIntPr<P, _64>(q, false) * tinv;
		}
		return qpow(x, ModIntPr<P, _64>::mod() - 2);
	}
}