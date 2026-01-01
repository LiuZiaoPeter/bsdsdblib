#pragma once

#include <vector>

#include "../../basics.hpp"

namespace numtheo_n {
	template<i32 P> MIP<P> qpow_signed(MIP<P> x, i64 y) {
		y %= MIP<P>::mod() - 1;
		if (y < 0) {
			y += MIP<P>::mod() - 1;
		}
		return qpow(x, y);
	}
	template<i32 P> void MIP<P>::lin_inv_preproc(u32 N) {
		MIP<P>::inv_v.resize(N + 1);
		MIP<P>::inv_v[0] = 0, MIP<P>::inv_v[1] = 1;
		for (u32 i = 2; i <= N; ++i) {
			MIP<P>::inv_v[i] = (-static_cast<MIP<P>>(MIP<P>::mod() / i))
					* MIP<P>::inv_v[MIP<P>::mod() % i];
		}
	}
	template<i32 P> void MIP<P>::O1inv_preproc() {
		MIP<P>::O1inv_mode = true;
		MIP<P>::cbrtP_log2 = std::bit_width(static_cast<u32>(std::cbrt(MIP<P>::mod()) + 1));
		MIP<P>::cbrtP = 1 << MIP<P>::cbrtP_log2;
		MIP<P>::cbrtP2 = 1 << (MIP<P>::cbrtP_log2 << 1);
		if (MIP<P>::inv_v.size() <= MIP<P>::cbrtP2) {
			lin_inv_preproc(MIP<P>::cbrtP2);
		}
		MIP<P>::farey_v.assign(MIP<P>::cbrtP2 + 1, std::make_pair(0, 0));
		for (u32 p = 0; p <= MIP<P>::cbrtP; ++p) {
			for (u32 q = (p == 1 ? 1 : p + 1); q <= MIP<P>::cbrtP; ++q) {
				u32 cur = (p << (MIP<P>::cbrtP_log2 << 1)) / q;
				if (MIP<P>::farey_v[cur].second == 0) {
					MIP<P>::farey_v[cur] = std::make_pair(p, q);
				}
			}
		}
		MIP<P>::farey_prec = MIP<P>::farey_v;
		for (auto it = MIP<P>::farey_prec.begin(); it != MIP<P>::farey_prec.end(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
		MIP<P>::farey_succ = MIP<P>::farey_v;
		for (auto it = MIP<P>::farey_succ.rbegin(); it != MIP<P>::farey_succ.rend(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
	}
	template<i32 P> MIP<P> inv(MIP<P> x) {
		if (x.val == 0) {
			throw std::invalid_argument(
				"numtheo_n::template<i32 P> MIP<P> inv(void) : getting inv(0)"
			);
		}
		if (x.val < MIP<P>::inv_v.size()) {
			return MIP<P>::inv_v[x.val];
		}
		if (MIP<P>::O1inv_mode == true) {
			u32 cur = static_cast<u32>(
				(static_cast<u64>(x.val) << (MIP<P>::cbrtP_log2 << 1))
				/ MIP<P>::mod()
			);
			u32 q;
			i32 t = static_cast<u32>(
				static_cast<i64>(x.val) * MIP<P>::farey_prec[cur].second
				- static_cast<i64>(MIP<P>::mod()) * MIP<P>::farey_prec[cur].first
			);
			if (std::abs(t) <= MIP<P>::cbrtP2) {
				q = MIP<P>::farey_prec[cur].second;
			} else {
				q = MIP<P>::farey_succ[cur].second;
				t = static_cast<u32>(
					static_cast<i64>(x.val) * q
					- static_cast<i64>(MIP<P>::mod()) * MIP<P>::farey_succ[cur].first
				);
			}
			MIP<P> tinv = ((t >= 0) ? MIP<P>::inv_v[t] : -MIP<P>::inv_v[-t]);
			return static_cast<MIP<P>>(q) * tinv;
		}
		return qpow(x, MIP<P>::mod() - 2);
	}
}