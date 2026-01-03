#pragma once

#include "../../basics.hpp"

namespace numtheo_n {
	template<i64 P, bool _64> void MIP<P, _64>::dis_ln_preproc(MIP<P, _64> g) {
		typename MIP<P, _64>::val_t sqrtP = static_cast<typename MIP<P, _64>::val_t>(std::sqrt(mod())) + 2;
		if (mpf.size() <= sqrtP) {
			euler_sieve(sqrtP);
		}
		std::vector<MIP<P, _64>> prs;
		for (typename MIP<P, _64>::val_t i = 2; i <= sqrtP; ++i) {
			if (mpf[i] == i) {
				prs.emplace_back(i, false);
			}
		}
		auto prln = dis_logs(g, prs);
		lesqrt_ln.resize(sqrtP + 1);
		for (typename MIP<P, _64>::val_t i = 0; i < prs.size(); ++i) {
			if (prln[i].has_value() == false) {
				throw std::invalid_argument(
					"template<i32 P> MIP<P, _64>::static void dis_ln_preproc(MIP<P, _64>) "
					": g not primitive root"
				);
			}
			lesqrt_ln[prs[i].val] = prln[i].value();
		}
		lesqrt_ln[1] = 0;
		for (typename MIP<P, _64>::val_t i = 4; i <= sqrtP; ++i) {
			lesqrt_ln[i] = lesqrt_ln[mpf[i]] + lesqrt_ln[i / mpf[i]];
			if (lesqrt_ln[i] >= mod() - 1) {
				lesqrt_ln[i] -= mod() - 1;
			}
		}
	}
	template<i64 P, bool _64> std::conditional_t<_64, u64, u32> fast_dis_ln(MIP<P, _64> x) {
		if (x.val < MIP<P, _64>::lesqrt_ln.size()) {
			return MIP<P, _64>::lesqrt_ln[x.val];
		}
		typename MIP<P, _64>::val_t k = MIP<P, _64>::mod() / x.val, r = MIP<P, _64>::mod() % x.val;
		if (r <= x.val - r) {
			typename MIP<P, _64>::mul_t ret = static_cast<typename MIP<P, _64>::mul_t>((MIP<P, _64>::mod() - 1) >> 1) + fast_dis_ln(MIP<P, _64>(r, false)) + (MIP<P, _64>::mod() - 1) - MIP<P, _64>::lesqrt_ln[k];
			if (ret >= MIP<P, _64>::mod() - 1) {
				ret -= MIP<P, _64>::mod() - 1;
				if (ret >= MIP<P, _64>::mod() - 1) {
					ret -= MIP<P, _64>::mod() - 1;
				}
			}
			return static_cast<typename MIP<P, _64>::val_t>(ret);
		} else {
			typename MIP<P, _64>::val_t ret = fast_dis_ln(MIP<P, _64>(x.val - r, false)) + (MIP<P, _64>::mod() - 1) - MIP<P, _64>::lesqrt_ln[k + 1];
			if (ret >= MIP<P, _64>::mod() - 1) {
				ret -= MIP<P, _64>::mod() - 1;
			}
			return ret;
		}
	}
}