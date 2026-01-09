// tested by: lg_P4777

#pragma once

#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../basics.hpp"

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
		using sT = i128::make_signed_t<T>;
		using uT = i128::up_t<T>;
		using usT = i128::up_t<sT>;
		if (a1 > a2) {
			std::swap(a1, a2), std::swap(p1, p2);
		}
		if ((a2 - a1) % std::gcd(p1, p2) != 0) {
			return std::nullopt;
		}
		auto [k1, k2] = exgcd(p1, p2);
		T g = std::gcd(p1, p2), l = std::lcm(p1, p2);
		sT g_s = static_cast<sT>(g), l_s = static_cast<sT>(l), mul = (a2 - a1) / g_s;
		k1 = static_cast<sT>(static_cast<usT>(k1) * mul % l_s);
		if (k1 < 0) {
			k1 += l;
		}
		k2 = static_cast<sT>(static_cast<usT>(-k2) * mul % l_s);
		if (k2 < 0) {
			k2 += l;
		}
		return std::make_pair((static_cast<uT>(k1) * p1 + a1) % l, l);
	}
	template<i128::liftable_unsigned T> std::optional<std::pair<T, T>> excrt(std::vector<T> a, std::vector<T> p) {
		if (a.size() != p.size() || a.empty() || p.empty()) {
			throw std::invalid_argument(
				static_cast<std::string>(__func__) + " : argument illegal"
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