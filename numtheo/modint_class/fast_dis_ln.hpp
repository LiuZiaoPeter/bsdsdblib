#pragma once

#include "../../basics.hpp"

namespace numtheo_n {
	template<i64 P, bool _64> void ModIntPr<P, _64>::dis_ln_preproc(ModIntPr<P, _64> g) {
		typename ModIntPr<P, _64>::val_t sqrtP = static_cast<typename ModIntPr<P, _64>::val_t>(std::sqrt(mod())) + 2;
		if (mpf.size() <= sqrtP) {
			euler_sieve(sqrtP);
		}
		std::vector<ModIntPr<P, _64>> prs;
		for (typename ModIntPr<P, _64>::val_t i = 2; i <= sqrtP; ++i) {
			if (mpf[i] == i) {
				prs.emplace_back(i, false);
			}
		}
		auto prln = dis_logs(g, prs);
		lesqrt_ln.resize(sqrtP + 1);
		for (typename ModIntPr<P, _64>::val_t i = 0; i < prs.size(); ++i) {
			if (prln[i].has_value() == false) {
				throw std::invalid_argument(
					"template<i32 P> ModIntPr<P, _64>::static void dis_ln_preproc(ModIntPr<P, _64>) "
					": g not primitive root"
				);
			}
			lesqrt_ln[prs[i].val] = prln[i].value();
		}
		lesqrt_ln[1] = 0;
		for (typename ModIntPr<P, _64>::val_t i = 4; i <= sqrtP; ++i) {
			lesqrt_ln[i] = lesqrt_ln[mpf[i]] + lesqrt_ln[i / mpf[i]];
			if (lesqrt_ln[i] >= mod() - 1) {
				lesqrt_ln[i] -= mod() - 1;
			}
		}
	}
	template<i64 P, bool _64> std::conditional_t<_64, u64, u32> fast_dis_ln(ModIntPr<P, _64> x) {
		if (x.val < ModIntPr<P, _64>::lesqrt_ln.size()) {
			return ModIntPr<P, _64>::lesqrt_ln[x.val];
		}
		typename ModIntPr<P, _64>::val_t k = ModIntPr<P, _64>::mod() / x.val, r = ModIntPr<P, _64>::mod() % x.val;
		if (r <= x.val - r) {
			typename ModIntPr<P, _64>::mul_t ret = static_cast<typename ModIntPr<P, _64>::mul_t>((ModIntPr<P, _64>::mod() - 1) >> 1) + fast_dis_ln(ModIntPr<P, _64>(r, false)) + (ModIntPr<P, _64>::mod() - 1) - ModIntPr<P, _64>::lesqrt_ln[k];
			if (ret >= ModIntPr<P, _64>::mod() - 1) {
				ret -= ModIntPr<P, _64>::mod() - 1;
				if (ret >= ModIntPr<P, _64>::mod() - 1) {
					ret -= ModIntPr<P, _64>::mod() - 1;
				}
			}
			return static_cast<typename ModIntPr<P, _64>::val_t>(ret);
		} else {
			typename ModIntPr<P, _64>::val_t ret = fast_dis_ln(ModIntPr<P, _64>(x.val - r, false)) + (ModIntPr<P, _64>::mod() - 1) - ModIntPr<P, _64>::lesqrt_ln[k + 1];
			if (ret >= ModIntPr<P, _64>::mod() - 1) {
				ret -= ModIntPr<P, _64>::mod() - 1;
			}
			return ret;
		}
	}
}