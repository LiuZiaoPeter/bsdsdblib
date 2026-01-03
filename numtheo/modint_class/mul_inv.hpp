#pragma once

#include <vector>

#include "../../basics.hpp"
#include "../../int128.hpp"

namespace numtheo_n {
	template<i128::signed_integral T, i64 P, bool _64> ModIntPr<P, _64> qpow_signed(ModIntPr<P, _64> x, T y) {
		using MIP = ModIntPr<P, _64>;
		y %= MIP::mod() - 1;
		if (y < 0) {
			y += MIP::mod() - 1;
		}
		return qpow(x, y);
	}
	template<i64 P, bool _64> template<i128::unsigned_integral T> void ModIntPr<P, _64>::lin_inv_preproc(T N) {
		using MIP = ModIntPr<P, _64>;
		MIP::inv_v.resize(N + 1);
		MIP::inv_v[0] = 0, MIP::inv_v[1] = 1;
		for (u32 i = 2; i <= N; ++i) {
			MIP::inv_v[i] = (-MIP(MIP::mod() / i, false)) * MIP::inv_v[MIP::mod() % i];
		}
	}
	template<i64 P, bool _64> void ModIntPr<P, _64>::O1inv_preproc() {
		using MIP = ModIntPr<P, _64>;
		MIP::O1inv_mode = true;
		MIP::cbrtP_log2 = std::bit_width(static_cast<MIP::val_t>(std::cbrt(MIP::mod()) + 1));
		MIP::cbrtP = static_cast<MIP::val_t>(1) << MIP::cbrtP_log2;
		MIP::cbrtP2 = static_cast<MIP::val_t>(1) << (MIP::cbrtP_log2 << 1);
		if (MIP::inv_v.size() <= MIP::cbrtP2) {
			lin_inv_preproc(MIP::cbrtP2);
		}
		MIP::farey_v.assign(MIP::cbrtP2 + 1, std::make_pair(0, 0));
		for (MIP::val_t p = 0; p <= MIP::cbrtP; ++p) {
			for (MIP::val_t q = (p == 1 ? 1 : p + 1); q <= MIP::cbrtP; ++q) {
				MIP::val_t cur = (p << (MIP::cbrtP_log2 << 1)) / q;
				if (MIP::farey_v[cur].second == 0) {
					MIP::farey_v[cur] = std::make_pair(p, q);
				}
			}
		}
		MIP::farey_prec = MIP::farey_v;
		for (auto it = MIP::farey_prec.begin(); it != MIP::farey_prec.end(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
		MIP::farey_succ = MIP::farey_v;
		for (auto it = MIP::farey_succ.rbegin(); it != MIP::farey_succ.rend(); ++it) {
			if (it->second == 0) {
				*it = *std::prev(it);
			}
		}
	}
	template<i64 P, bool _64> ModIntPr<P, _64> inv(ModIntPr<P, _64> x) {
		using MIP = ModIntPr<P, _64>;
		using val_t = typename MIP::val_t;
		using sval_t = i128::make_signed_t<val_t>;
		using mul_t = typename MIP::mul_t;
		using smul_t = i128::make_signed_t<mul_t>;
		if (x.val == 0) {
			throw std::invalid_argument("numtheo_n::template<i64 P, bool _64> ModIntPr<P, _64> inv(void) : getting inv(0)");
		}
		if (x.val < MIP::inv_v.size()) {
			return MIP::inv_v[x.val];
		}
		if (MIP::O1inv_mode == true) {
			val_t cur = static_cast<val_t>((static_cast<mul_t>(x.val) << (MIP::cbrtP_log2 << 1)) / MIP::mod());
			val_t q;
			sval_t t = static_cast<sval_t>(static_cast<smul_t>(x.val) * MIP::farey_prec[cur].second - static_cast<smul_t>(MIP::mod()) * MIP::farey_prec[cur].first);
			if (static_cast<val_t>(std::abs(t)) <= MIP::cbrtP2) {
				q = MIP::farey_prec[cur].second;
			} else {
				q = MIP::farey_succ[cur].second;
				t = static_cast<sval_t>(static_cast<smul_t>(x.val) * q - static_cast<smul_t>(MIP::mod()) * MIP::farey_succ[cur].first);
			}
			MIP tinv = ((t >= 0) ? MIP::inv_v[t] : -MIP::inv_v[-t]);
			return MIP(q, false) * tinv;
		}
		return qpow(x, MIP::mod() - 2);
	}
}