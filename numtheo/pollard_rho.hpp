// tested by lg_P4718

#pragma once

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

#include "../basics.hpp"

namespace numtheo {
	template<i128::liftable_unsigned T> bool miller_rabin(T x) {
		if (!(x & 1)) {
			return x == 2;
		}
		T t = __builtin_ctzll(x - 1), u = (x - 1) >> t;
		struct modx {
			T val, _mod;
			modx &operator*=(modx _x) {
				val = static_cast<T>(static_cast<i128::up_t<T>>(val) * _x.val % _mod);
				return *this;
			}
		};
		std::vector<T> a_list;
		if (std::is_same_v<T, u64>) {
			a_list = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
		} else {
			a_list = {2, 7, 61};
		}
		for (T a : a_list) {
			a %= x;
			if (a == 0) {
				continue;
			}
			if (std::gcd(a, x) != 1) {
				return false;
			}
			modx v = qpow(modx{a, x}, u, modx{1, x});
			if (v.val == 1) {
				continue;
			}
			bool nfound = false;
			for (u32 s = 0; s < t; ++s) {
				if (v.val == x - 1) {
					nfound = true;
					break;
				}
				v *= v;
			}
			if (nfound == false) {
				return false;
			}
		}
		return true;
	}
	template<i128::liftable_unsigned T> T pollard_rho(T x) {
		if (!(x & 1)) {
			return 2;
		}
		T c;
		std::conditional_t<std::is_same_v<T, u64>, std::mt19937_64, std::mt19937> rndc(std::random_device{}());
		while (true) {
			c = static_cast<T>(rndc()) % (x - 1) + 1;
			auto f = [c, x](T _x)->T {
				return static_cast<T>((static_cast<i128::up_t<T>>(_x) * _x + c) % x);
			};
			T s = 0, t = 0, prod = 1;
			for (u32 k = 1;; ++k) {
				u32 counter = 0;
				for (u32 step = 0; !(step >> k); ++step) {
					t = f(t);
					prod = static_cast<T>(
						(static_cast<i128::up_t<T>>(prod) * (s < t ? t - s : s - t)) % x
					);
					++counter;
					if (counter == 128) {
						T g = std::gcd(prod, x);
						if (g == x) {
							goto newc;
						}
						if (g != 1) {
							return g;
						}
						counter = 0;
					}
				}
				T g = std::gcd(prod, x);
				if (g == x) {
					goto newc;
				}
				if (g != 1) {
					return g;
				}
				s = t;
			}
			newc:;
		}
	}
	template<i128::liftable_unsigned T> void prime_factors(T x, std::vector<T> &ret) {
		if (x <= 1) {
			return;
		}
		if (miller_rabin(x)) {
			ret.emplace_back(x);
			return;
		}
		T fact = pollard_rho(x);
		prime_factors(fact, ret);
		prime_factors(x / fact, ret);
		return;
	}
	template<i128::liftable_unsigned T> std::vector<std::pair<T, u32>> prime_factors(T x) {
		std::vector<T> ret1;
		prime_factors(x, ret1);
		std::sort(ret1.begin(), ret1.end());
		std::vector<std::pair<T, u32>> ret;
		for (T i : ret1) {
			if (ret.empty() || i != ret.back().first) {
				ret.emplace_back(i, 1);
			} else {
				++ret.back().second;
			}
		}
		return ret;
	}
}