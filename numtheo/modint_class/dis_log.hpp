#pragma once

#include <cmath>
#include <optional>
#include <unordered_map>
#include <vector>

#include "../../basics.hpp"
#include "../prod_funcs.hpp"

namespace numtheo_n {
	template<i64 P, bool _64> std::optional<std::conditional_t<_64, u64, u32>> dis_log(MIP<P, _64> a, MIP<P, _64> b) {
		typename MIP<P, _64>::val_t B = static_cast<typename MIP<P, _64>::val_t>(std::sqrt(MIP<P, _64>::mod())) + 2;
		MIP<P, _64> a_to_y = 1;
		std::unordered_map<typename MIP<P, _64>::val_t, typename MIP<P, _64>::val_t> bay2y;
		for (typename MIP<P, _64>::val_t y = 0; y < B; ++y, a_to_y *= a) {
			bay2y[(b * a_to_y).val] = y;
		}
		MIP<P, _64> a_to_B = a_to_y, a_to_B_to_x = a_to_B;
		for (typename MIP<P, _64>::val_t x = 1; x <= B; ++x, a_to_B_to_x *= a_to_B) {
			if (bay2y.find(a_to_B_to_x.val) != bay2y.end()) {
				return static_cast<typename MIP<P, _64>::mul_t>(B) * x - bay2y[a_to_B_to_x.val];
			}
		}
		return std::nullopt;
	}
	template<i64 P, bool _64> std::vector<std::optional<std::conditional_t<_64, u64, u32>>> dis_logs(MIP<P, _64> a, std::vector<MIP<P, _64>> b) {
		typename MIP<P, _64>::val_t B = static_cast<typename MIP<P, _64>::val_t>(std::sqrt(MIP<P, _64>::mod() / b.size())) + 2;
		typename MIP<P, _64>::val_t xlim = MIP<P, _64>::mod() / B + 3;
		MIP<P, _64> a_to_B = qpow(a, B), a_to_B_to_x = a_to_B;
		std::unordered_map<typename MIP<P, _64>::val_t, typename MIP<P, _64>::val_t> aBx2x;
		for (typename MIP<P, _64>::val_t x = 1; x <= xlim; ++x, a_to_B_to_x *= a_to_B) {
			if (aBx2x.find(a_to_B_to_x.val) == aBx2x.end()) {
				aBx2x[a_to_B_to_x.val] = x;
			}
		}
		std::vector<std::optional<typename MIP<P, _64>::val_t>> ret(b.size(), std::nullopt);
		MIP<P, _64> a_to_y = 1;
		for (typename MIP<P, _64>::val_t y = 0; y < B; ++y, a_to_y *= a) {
			for (typename MIP<P, _64>::val_t i = 0; i < ret.size(); ++i) {
				typename MIP<P, _64>::val_t bayv = (b[i] * a_to_y).val;
				if (aBx2x.find(bayv) == aBx2x.end()) {
					continue;
				}
				typename MIP<P, _64>::val_t cura = static_cast<typename MIP<P, _64>::val_t>(static_cast<typename MIP<P, _64>::mul_t>(aBx2x[bayv]) * B - y);
				if (ret[i].has_value() == false) {
					ret[i] = cura;
				} else {
					ret[i] = std::min(ret[i].value(), cura);
				}
			}
		}
		return ret;
	}
	template<i64 P, bool _64> std::conditional_t<_64, u64, u32> ord(MIP<P, _64> x) {
		return dis_log(x, MIP<P, _64>(1, false));
	}
	template<i32 P> std::optional<u32> dis_log(MI<P> a, MI<P> b) {
		u32 B = static_cast<u32>(std::sqrt(MI<P>::mod())) + 1;
		MI<P> a_to_B = qpow(a, B), a_to_Bx = a_to_B;
		std::unordered_map<u32, std::pair<u32, u32>> aBx2x;
		for (u32 x = 1; x <= B; ++x, a_to_Bx *= a_to_B) {
			if (aBx2x[a_to_Bx.val].first == 0) {
				aBx2x[a_to_Bx.val].first = x;
			} else if (aBx2x[a_to_Bx.val].second == 0) {
				aBx2x[a_to_Bx.val].second = x;
			}
		}
		MI<P> a_to_y = 1;
		u32 ret = -1;
		for (u32 y = 0; y < B; ++y, a_to_y *= a) {
			u32 curv = (b * a_to_y).val;
			if (aBx2x.find(curv) == aBx2x.end()) {
				continue;
			}
			for (u32 x : {aBx2x[curv].first, aBx2x[curv].second}) {
				if (x == 0) {
					continue;
				}
				u32 cur = static_cast<u32>(static_cast<u64>(B) * x - y);
				if (qpow(a, cur) == b) {
					ret = std::min(ret, cur);
					break;
				}
			}
		}
		if (ret == static_cast<u32>(-1)) {
			return std::nullopt;
		}
		return ret;
	}
	template<i32 P> std::vector<std::optional<u32>> dis_logs(MI<P> a, std::vector<MI<P>> b) {
		u32 B = static_cast<u32>(std::sqrt(phi(MI<P>::mod()) / b.size())) + 2;
		u32 xlim = MI<P>::mod() / B + 3;
		MI<P> a_to_B = qpow(a, B), a_to_Bx = a_to_B;
		std::unordered_map<u32, std::pair<u32, u32>> aBx2x;
		for (u32 x = 1; x <= xlim; ++x, a_to_Bx *= a_to_B) {
			if (aBx2x[a_to_Bx.val].first == 0) {
				aBx2x[a_to_Bx.val].first = x;
			} else if (aBx2x[a_to_Bx.val].second == 0) {
				aBx2x[a_to_Bx.val].second = x;
			}
		}
		MI<P> a_to_y = 1;
		std::vector<std::optional<u32>> ret(b.size(), std::nullopt);
		for (u32 y = 0; y < B; ++y, a_to_y *= a) {
			for (u32 i = 0; i < B; ++i) {
				u32 curv = (b[i] * a_to_y).val;
				if (aBx2x.find(curv) == aBx2x.end()) {
					continue;
				}
				for (u32 x : {aBx2x[curv].first, aBx2x[curv].second}) {
					if (x == 0) {
						continue;
					}
					u32 cur = static_cast<u32>(static_cast<u64>(B) * x - y);
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
	template<i32 P> std::optional<u32> ord(MI<P> x) {
		if (std::gcd(x.val, MI<P>::mod()) != 1) {
			return std::nullopt;
		}
		MIP<-1073741824>::set_mod(MI<P>::mod());
		return dis_log(MIP<-1073741824>(x.val, false), MIP<-1073741824>(1, false));
	}
}