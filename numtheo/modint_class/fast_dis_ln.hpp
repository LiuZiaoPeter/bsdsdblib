#pragma once

#include "../../basics.hpp"

namespace numtheo_n {
	template<i64 P, bool _64> void ModIntPr<P, _64>::dis_ln_preproc(ModIntPr<P, _64> g) {
		using MIP = ModIntPr<P, _64>;
		using val_t = typename MIP::val_t;
		typename ModIntPr<P, _64>::val_t sqrtP = static_cast<typename ModIntPr<P, _64>::val_t>(std::sqrt(mod())) + 2;
		if (mpf.size() <= sqrtP) {
			euler_sieve(sqrtP);
		}
		std::vector<MIP> prs;
		for (val_t i = 2; i <= sqrtP; ++i) {
			if (mpf[i] == i) {
				prs.emplace_back(i, false);
			}
		}
		auto prln = dis_logs(g, prs);
		lesqrt_ln.resize(sqrtP + 1);
		for (val_t i = 0; i < prs.size(); ++i) {
			if (prln[i].has_value() == false) {
				throw std::invalid_argument(
					"template<i64 P, bool _64> ModIntPr<P, _64>::static void dis_ln_preproc"
					"(ModIntPr<P, _64>) : g not primitive root"
				);
			}
			lesqrt_ln[prs[i].val] = prln[i].value();
		}
		lesqrt_ln[1] = 0;
		for (val_t i = 4; i <= sqrtP; ++i) {
			lesqrt_ln[i] = lesqrt_ln[mpf[i]] + lesqrt_ln[i / mpf[i]];
			if (lesqrt_ln[i] >= mod() - 1) {
				lesqrt_ln[i] -= mod() - 1;
			}
		}
	}
	template<i64 P, bool _64> std::conditional_t<_64, u64, u32> fast_dis_ln(ModIntPr<P, _64> x) {
		using MIP = ModIntPr<P, _64>;
		using val_t = typename MIP::val_t;
		using mul_t = typename MIP::mul_t;
		if (x.val < MIP::lesqrt_ln.size()) {
			return MIP::lesqrt_ln[x.val];
		}
		val_t k = MIP::mod() / x.val, r = MIP::mod() % x.val;
		if (r <= x.val - r) {
			mul_t ret = static_cast<mul_t>((MIP::mod() - 1) >> 1) + fast_dis_ln(MIP(r, false))
					+ (MIP::mod() - 1) - MIP::lesqrt_ln[k];
			if (ret >= MIP::mod() - 1) {
				ret -= MIP::mod() - 1;
				if (ret >= MIP::mod() - 1) {
					ret -= MIP::mod() - 1;
				}
			}
			return static_cast<val_t>(ret);
		} else {
			val_t ret = fast_dis_ln(MIP(x.val - r, false)) + (MIP::mod() - 1) - MIP::lesqrt_ln[k + 1];
			if (ret >= MIP::mod() - 1) {
				ret -= MIP::mod() - 1;
			}
			return ret;
		}
	}
}