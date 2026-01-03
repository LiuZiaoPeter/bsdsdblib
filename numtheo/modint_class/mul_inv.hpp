#pragma once

#include <vector>

#include "../../basics.hpp"
#include "../../int128.hpp"

namespace numtheo_n {
	template<i128::signed_integral T, i64 P, bool _64> MIP<P, _64> qpow_signed(MIP<P, _64> x, T y) {
		y %= MIP<P, _64>::mod() - 1;
		if (y < 0) {
			y += MIP<P, _64>::mod() - 1;
		}
		return qpow(x, y);
	}
	template<i64 P, bool _64> template<i128::unsigned_integral T> void MIP<P, _64>::lin_inv_preproc(T N) {
		MIP<P, _64>::inv_v.resize(N + 1);
		MIP<P, _64>::inv_v[0] = 0, MIP<P, _64>::inv_v[1] = 1;
		for (u32 i = 2; i <= N; ++i) {
			MIP<P, _64>::inv_v[i] = (-MIP<P, _64>(MIP<P, _64>::mod() / i, false)) * MIP<P, _64>::inv_v[MIP<P, _64>::mod() % i];
		}
	}
	template<i64 P, bool _64> void MIP<P, _64>::O1inv_preproc() {
		MIP<P, _64>::O1inv_mode = true;
		MIP<P, _64>::cbrtP_log2 = std::bit_width(
			static_cast<MIP<P, _64>::val_t>(std::cbrt(MIP<P, _64>::mod()) + 1)
		);
		MIP<P, _64>::cbrtP = static_cast<MIP<P, _64>::val_t>(1) << MIP<P, _64>::cbrtP_log2;
		MIP<P, _64>::cbrtP2 = static_cast<MIP<P, _64>::val_t>(1) << (MIP<P, _64>::cbrtP_log2 << 1);
		if (MIP<P, _64>::inv_v.size() <= MIP<P, _64>::cbrtP2) {
			lin_inv_preproc(MIP<P, _64>::cbrtP2);
		}
		MIP<P, _64>::farey_v.assign(MIP<P, _64>::cbrtP2 + 1, std::make_pair(0, 0));
		for (MIP<P, _64>::val_t p = 0; p <= MIP<P, _64>::cbrtP; ++p) {
			for (MIP<P, _64>::val_t q = (p == 1 ? 1 : p + 1); q <= MIP<P, _64>::cbrtP; ++q) {
				MIP<P, _64>::val_t cur = (p << (MIP<P, _64>::cbrtP_log2 << 1)) / q;
				if (MIP<P, _64>::farey_v[cur].second == 0) {
					MIP<P, _64>::farey_v[cur] = std::make_pair(p, q);
				}
			}
		}
		MIP<P, _64>::farey_prec = MIP<P, _64>::farey_v;
		for (auto it = MIP<P, _64>::farey_prec.begin(); it != MIP<P, _64>::farey_prec.end(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
		MIP<P, _64>::farey_succ = MIP<P, _64>::farey_v;
		for (auto it = MIP<P, _64>::farey_succ.rbegin(); it != MIP<P, _64>::farey_succ.rend(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
	}
	template<i64 P, bool _64> MIP<P, _64> inv(MIP<P, _64> x) {
		if (x.val == 0) {
			throw std::invalid_argument(
				"numtheo_n::template<i32 P> MIP<P, _64> inv(void) : getting inv(0)"
			);
		}
		if (x.val < MIP<P, _64>::inv_v.size()) {
			return MIP<P, _64>::inv_v[x.val];
		}
		if (MIP<P, _64>::O1inv_mode == true) {
			typename MIP<P, _64>::val_t cur = static_cast<MIP<P, _64>::val_t>((static_cast<MIP<P, _64>::mul_t>(x.val) << (MIP<P, _64>::cbrtP_log2 << 1)) / MIP<P, _64>::mod());
			typename MIP<P, _64>::val_t q;
			i128::make_signed_t<typename MIP<P, _64>::val_t> t = static_cast<i128::make_signed_t<typename MIP<P, _64>::val_t>>(static_cast<i128::make_signed_t<typename MIP<P, _64>::mul_t>>(x.val) * MIP<P, _64>::farey_prec[cur].second - static_cast<i128::make_signed_t<typename MIP<P, _64>::mul_t>>(MIP<P, _64>::mod()) * MIP<P, _64>::farey_prec[cur].first);
			if (static_cast<typename MIP<P, _64>::val_t>(std::abs(t)) <= MIP<P, _64>::cbrtP2) {
				q = MIP<P, _64>::farey_prec[cur].second;
			} else {
				q = MIP<P, _64>::farey_succ[cur].second;
				t = static_cast<i128::make_signed_t<typename MIP<P, _64>::val_t>>(static_cast<i128::make_signed_t<typename MIP<P, _64>::mul_t>>(x.val) * q - static_cast<i128::make_signed_t<typename MIP<P, _64>::mul_t>>(MIP<P, _64>::mod()) * MIP<P, _64>::farey_succ[cur].first);
			}
			MIP<P, _64> tinv = ((t >= 0) ? MIP<P, _64>::inv_v[t] : -MIP<P, _64>::inv_v[-t]);
			return MIP<P, _64>(q, false) * tinv;
		}
		return qpow(x, MIP<P, _64>::mod() - 2);
	}
}