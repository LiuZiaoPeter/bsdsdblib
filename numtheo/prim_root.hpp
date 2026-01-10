#pragma once

#include <numeric>
#include <optional>
#include <vector>

#include "../basics.hpp"
#include "pollard_rho.hpp"
#include "prod_funcs.hpp"

namespace numtheo {
	template<i128::liftable_unsigned T> bool is_prim_root(T x, T m, T phim, const std::vector<std::pair<T, u32>> &phim_fact) {
		using MI = ModInt<modint_inner + 5, std::is_same_v<T, u64>>;
		MI::set_mod(m);
		if (std::gcd(x, m) != 1) {
			return false;
		}
		for (auto i : phim_fact) {
			if (qpow(MI(x, false), phim / i.first, MI(1, false)).value() == 1) {
				return false;
			}
		}
		return true;
	}
	template<i128::liftable_unsigned T> bool is_prim_root(T x, T m) {
		T phim = phi(m);
		auto phim_fact = prime_factors(phim);
		return is_prim_root(x, m, phim, phim_fact);
	}
	template<i128::liftable_unsigned T> std::optional<T> min_prim_root(T m, T phim, const std::vector<std::pair<T, u32>> &phim_fact) {
		if (m == 2) {
			return 1;
		}
		if (m == 4) {
			return 3;
		}
		// check
		T _m = m;
		if (!(_m & 1)) {
			_m /= 2;
		}
		if (!(_m & 1) || prime_factors(_m).size() != 1) {
			return std::nullopt;
		}
		// calculate
		for (T i = 1;; ++i) {
			if (is_prim_root(i, m, phim, phim_fact)) {
				return i;
			}
		}
	}
	template<i128::liftable_unsigned T> std::optional<T> min_prim_root(T m) {
		return min_prim_root(m, phi(m), prime_factors(phi(m)));
	}
	template<i128::liftable_unsigned T> std::vector<T> prim_roots(T m) {
		if (m == 2) {
			return {1};
		}
		if (m == 4) {
			return {3};
		}
		T phim = phi(m);
		auto phim_fact = prime_factors(phim);
		std::optional<T> min_pr = min_prim_root(m, phim, phim_fact);
		if (min_pr.has_value() == false) {
			return std::vector<T>();
		}
		std::vector<T> ret;
		T g_to_k = 1;
		for (T k = 1; k <= phim; ++k) {
			g_to_k = g_to_k * min_pr.value() % m;
			if (std::gcd(k, phim) == 1) {
				ret.emplace_back(g_to_k);
			}
		}
		std::sort(ret.begin(), ret.end());
		return ret;
	}
}