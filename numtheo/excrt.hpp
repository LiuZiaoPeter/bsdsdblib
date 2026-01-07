// tested by: lg_P4777

#pragma once

#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../basics.hpp"

namespace numtheo_n {
        std::pair<i64, i64> exgcd(u64 a, u64 b) {
		if (b == 0) {
			return std::make_pair(1, 0);
		}
		u64 k = a / b;
		std::pair<i64, i64> pr = exgcd(b, a % b);
		return std::make_pair(pr.second, pr.first - k * pr.second);
	}
	std::optional<std::pair<u64, u64>> excrt(u64 a1, u64 p1, u64 a2, u64 p2) {
		if (a1 > a2) {
			std::swap(a1, a2), std::swap(p1, p2);
		}
		if ((a2 - a1) % std::gcd(p1, p2) != 0) {
			return std::nullopt;
		}
		auto [k1, k2] = exgcd(p1, p2);
		u64 g = std::gcd(p1, p2), l = std::lcm(p1, p2);
		i64 g_s = g, l_s = l, mul = (a2 - a1) / g_s;
		k1 = static_cast<i64>(static_cast<__int128_t>(k1) * mul % l_s);
		if (k1 < 0) {
			k1 += l;
		}
		k2 = static_cast<i64>(static_cast<__int128_t>(-k2) * mul % l_s);
		if (k2 < 0) {
			k2 += l;
		}
		// assert((__uint128_t(k1) * p1 + a1) % l == (__uint128_t(k2) * p2 + a2) % l);
		return std::make_pair((static_cast<__uint128_t>(k1) * p1 + a1) % l, l);
	}
	std::optional<std::pair<u64, u64>> excrt(std::vector<u64> a, std::vector<u64> p) {
		if (a.size() != p.size() || a.empty() || p.empty()) {
			throw std::invalid_argument(
				"numtheo_n::std::optional<std::pair<u64, u64>> excrt"
				"(std::vector<u64>, std::vector<u64>) : argument illegal"
			);
		}
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