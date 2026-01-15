#pragma once

#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../basics.hpp"
#include "modint.hpp"

namespace numtheo {
        template<i128::unsigned_integral T> std::pair<i128::make_signed_t<T>, i128::make_signed_t<T>> exgcd(T a, T b) {
		using sT = i128::make_signed_t<T>;
		if (b == 0) {
			return std::make_pair(1, 0);
		}
		u64 k = a / b;
		std::pair<sT, sT> pr = exgcd(b, a % b);
		return std::make_pair(pr.second, pr.first - k * pr.second);
	}
	template<i128::liftable_unsigned T> std::optional<std::pair<T, T>> excrt(T a1, T p1, T a2, T p2) {
		using MI = ModInt<modint_inner + 2, std::is_same_v<T, u64>>; // mod : lcm(p1, p2)
		if (a1 > a2) {
			std::swap(a1, a2), std::swap(p1, p2);
		}
		if ((a2 - a1) % std::gcd(p1, p2) != 0) {
			return std::nullopt;
		}
		auto k1 = exgcd(p1, p2).first;
		T g = std::gcd(p1, p2), l = std::lcm(p1, p2), mul = (a2 - a1) / g;
		MI::set_mod(l);
		MI k1m = MI(k1) * MI(mul, false);
		return std::make_pair((k1m * MI(p1, false) + MI(a1, false)).value(), l);
	}
	template<i128::liftable_unsigned T> std::optional<std::pair<T, T>> excrt(std::vector<T> a, std::vector<T> p) {
		assure(a.size() == p.size(), "size(a)={} while size(p)={}", a.size(), p.size());
		assure(!a.empty(), "a is empty");
		assure(!p.empty(), "p is empty");
		for (auto ita = a.begin(), itp = p.begin(); ita != prev(a.end()); ++ita, ++itp) {
			auto res = excrt(*ita, *itp, *next(ita), *next(itp));
			if (res.has_value() == false) {
				return std::nullopt;
			}
			std::tie(*next(ita), *next(itp)) = res.value();
		}
		return std::make_pair(a.back(), p.back());
	}
}
