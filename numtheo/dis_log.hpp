// tested by : lg_P3846, lg_P4195, lg_U648171, lg_U648236

#pragma once

#include <cmath>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <vector>

#include "../basics.hpp"
#include "modint.hpp"
#include "prod_funcs.hpp"

namespace numtheo {
	template<i64 P, bool _64>
	std::optional<std::conditional_t<_64, u64, u32>> dis_log(ModIntPr<P, _64> a, ModIntPr<P, _64> b) {
		using MIP = ModIntPr<P, _64>;
		using val_t = MIP::val_t;
		using mul_t = MIP::mul_t;
		u32 B = static_cast<u32>(std::sqrt(MIP::mod())) + 2;
		MIP a_to_y = 1;
		std::unordered_map<val_t, val_t> bay2y;
		for (u32 y = 0; y < B; ++y, a_to_y *= a) {
			bay2y[(b * a_to_y).value()] = y;
		}
		MIP a_to_B = a_to_y, a_to_B_to_x = a_to_B;
		for (u32 x = 1; x <= B; ++x, a_to_B_to_x *= a_to_B) {
			if (bay2y.find(a_to_B_to_x.value()) != bay2y.end()) {
				return static_cast<mul_t>(B) * x - bay2y[a_to_B_to_x.value()];
			}
		}
		return std::nullopt;
	}
	template<i64 P, bool _64> std::vector<std::optional<std::conditional_t<_64, u64, u32>>>
	dis_logs(ModIntPr<P, _64> a, const std::vector<ModIntPr<P, _64>> &b) {
		using MIP = ModIntPr<P, _64>;
		using val_t = MIP::val_t;
		using mul_t = MIP::mul_t;
		u32 B = static_cast<u32>(std::sqrt(MIP::mod() / b.size())) + 2;
		u32 xlim = static_cast<u32>(MIP::mod() / B) + 3;
		MIP a_to_B = qpow(a, B), a_to_B_to_x = a_to_B;
		std::unordered_map<val_t, val_t> aBx2x;
		for (u32 x = 1; x <= xlim; ++x, a_to_B_to_x *= a_to_B) {
			if (aBx2x.find(a_to_B_to_x.value()) == aBx2x.end()) {
				aBx2x[a_to_B_to_x.value()] = x;
			}
		}
		std::vector<std::optional<val_t>> ret(b.size(), std::nullopt);
		MIP a_to_y = 1;
		for (u32 y = 0; y < B; ++y, a_to_y *= a) {
			for (u32 i = 0; i < ret.size(); ++i) {
				val_t bayv = (b[i] * a_to_y).value();
				if (aBx2x.find(bayv) == aBx2x.end()) {
					continue;
				}
				val_t cura = static_cast<val_t>(static_cast<mul_t>(aBx2x[bayv]) * B - y);
				if (ret[i].has_value() == false) {
					ret[i] = cura;
				} else {
					ret[i] = std::min(ret[i].value(), cura);
				}
			}
		}
		return ret;
	}
	template<i64 P, bool _64> std::conditional_t<_64, u64, u32> ord(ModIntPr<P, _64> x) {
		return dis_log(x, ModIntPr<P, _64>(1, false));
	}

	template<i64 P, bool _64>
	std::optional<std::conditional_t<_64, u64, u32>> dis_log(ModInt<P, _64> a, ModInt<P, _64> b) {
		using MI = ModInt<P, _64>;
		using val_t = MI::val_t;
		using mul_t = MI::mul_t;
		u32 B = static_cast<u32>(std::sqrt(MI::mod())) + 1;
		MI a_to_B = qpow(a, B), a_to_Bx = a_to_B;
		std::unordered_map<val_t, std::pair<val_t, val_t>> aBx2x;
		for (u32 x = 1; x <= B; ++x, a_to_Bx *= a_to_B) {
			if (aBx2x[a_to_Bx.value()].first == 0) {
				aBx2x[a_to_Bx.value()].first = x;
			} else if (aBx2x[a_to_Bx.value()].second == 0) {
				aBx2x[a_to_Bx.value()].second = x;
			}
		}
		MI a_to_y = 1;
		std::optional<val_t> ret = std::nullopt;
		for (u32 y = 0; y < B; ++y, a_to_y *= a) {
			val_t curv = (b * a_to_y).value();
			if (aBx2x.find(curv) == aBx2x.end()) {
				continue;
			}
			for (val_t x : {aBx2x[curv].first, aBx2x[curv].second}) {
				if (x == 0) {
					continue;
				}
				val_t cur = static_cast<val_t>(static_cast<mul_t>(B) * x - y);
				if (qpow(a, cur) == b) {
					if (ret.has_value()) {
						ret = std::min(ret.value(), cur);
					} else {
						ret = cur;
					}
					break;
				}
			}
		}
		return ret;
	}
	template<i64 P, bool _64> std::vector<std::optional<std::conditional_t<_64, u64, u32>>>
	dis_logs(ModInt<P, _64> a, const std::vector<ModInt<P, _64>> &b) {
		using MI = ModInt<P, _64>;
		using val_t = MI::val_t;
		using mul_t = MI::mul_t;
		u32 B = static_cast<u32>(std::sqrt(phi(MI::mod()) / b.size())) + 2;
		u32 xlim = static_cast<u32>(MI::mod() / B) + 3;
		MI a_to_B = qpow(a, B), a_to_Bx = a_to_B;
		std::unordered_map<val_t, std::pair<val_t, val_t>> aBx2x;
		for (u32 x = 1; x <= xlim; ++x, a_to_Bx *= a_to_B) {
			if (aBx2x[a_to_Bx.value()].first == 0) {
				aBx2x[a_to_Bx.value()].first = x;
			} else if (aBx2x[a_to_Bx.value()].second == 0) {
				aBx2x[a_to_Bx.value()].second = x;
			}
		}
		MI a_to_y = 1;
		std::vector<std::optional<val_t>> ret(b.size(), std::nullopt);
		for (u32 y = 0; y < B; ++y, a_to_y *= a) {
			for (u32 i = 0; i < ret.size(); ++i) {
				val_t curv = (b[i] * a_to_y).value();
				if (aBx2x.find(curv) == aBx2x.end()) {
					continue;
				}
				for (val_t x : {aBx2x[curv].first, aBx2x[curv].second}) {
					if (x == 0) {
						continue;
					}
					val_t cur = static_cast<val_t>(static_cast<mul_t>(B) * x - y);
					if (qpow(a, cur) == b[i]) {
						if (ret[i].has_value() == false) {
							ret[i] = cur;
						} else {
							ret[i] = std::min(ret[i].value(), cur);
						}
						break;
					}
				}
			}
		}
		return ret;
	}
	template<i64 P, bool _64> std::optional<std::conditional_t<_64, u64, u32>> ord(ModInt<P, _64> x) {
		using MI = ModInt<P, _64>;
		if (std::gcd(x.value(), MI::mod()) != 1) {
			return std::nullopt;
		}
		using tmpMIP = ModIntPr<-1073741824, _64>;
		tmpMIP::set_mod(MI::mod());
		return dis_log(tmpMIP(x.value(), false), tmpMIP(1, false));
	}
}