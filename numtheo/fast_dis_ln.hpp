#pragma once

#include <vector>

#include "../basics.hpp"
#include "modint.hpp"
#include "dis_log.hpp"
#include "euler_sieve.hpp"

namespace fast_dis_ln_hpp {
	template<i64 P, bool _64> struct aux {
		using val_t = std::conditional_t<_64, u64, u32>;
		inline static std::vector<val_t> lesqrt_ln;
	};
}

namespace numtheo {
	template<i64 P, bool _64> void dis_ln_preproc(ModIntPr<P, _64> g) {
		using MIP = ModIntPr<P, _64>;
		using aux = fast_dis_ln_hpp::aux<P, _64>;
		u32 sqrtP = static_cast<u32>(std::sqrt(MIP::mod())) + 2;
		if (mpf.size() <= sqrtP) {
			euler_sieve(sqrtP);
		}
		std::vector<MIP> prs;
		for (u32 i = 2; i <= sqrtP; ++i) {
			if (mpf[i] == i) {
				prs.emplace_back(i, false);
			}
		}
		auto prln = dis_logs(g, prs);
		aux::lesqrt_ln.resize(sqrtP + 1);
		for (u32 i = 0; i < prs.size(); ++i) {
			if (prln[i].has_value() == false) {
				throw std::invalid_argument(
					__func_str__ + " : g not primitive root"
				);
			}
			aux::lesqrt_ln[prs[i].value()] = prln[i].value();
		}
		aux::lesqrt_ln[1] = 0;
		for (u32 i = 4; i <= sqrtP; ++i) {
			aux::lesqrt_ln[i] = aux::lesqrt_ln[mpf[i]] + aux::lesqrt_ln[i / mpf[i]];
			if (aux::lesqrt_ln[i] >= MIP::mod() - 1) {
				aux::lesqrt_ln[i] -= MIP::mod() - 1;
			}
		}
	}
	template<i64 P, bool _64> std::conditional_t<_64, u64, u32> fast_dis_ln(ModIntPr<P, _64> x) {
		using MIP = ModIntPr<P, _64>;
		using aux = fast_dis_ln_hpp::aux<P, _64>;
		using val_t = MIP::val_t;
		using mul_t = MIP::mul_t;
		if (x.value() < aux::lesqrt_ln.size()) {
			return aux::lesqrt_ln[x.value()];
		}
		val_t k = MIP::mod() / x.value(), r = MIP::mod() % x.value();
		if (r <= x.value() - r) {
			mul_t ret = static_cast<mul_t>((MIP::mod() - 1) >> 1) + fast_dis_ln(MIP(r, false))
					+ (MIP::mod() - 1) - aux::lesqrt_ln[k];
			if (ret >= MIP::mod() - 1) {
				ret -= MIP::mod() - 1;
				if (ret >= MIP::mod() - 1) {
					ret -= MIP::mod() - 1;
				}
			}
			return static_cast<val_t>(ret);
		} else {
			val_t ret = fast_dis_ln(MIP(x.value() - r, false)) + (MIP::mod() - 1) - aux::lesqrt_ln[k + 1];
			if (ret >= MIP::mod() - 1) {
				ret -= MIP::mod() - 1;
			}
			return ret;
		}
	}
}
