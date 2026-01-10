#pragma once

#include <utility>
#include <vector>

#include "../basics.hpp"

namespace numtheo {
	template<i128::unsigned_integral T> std::vector<T> cont_frac(T a, T b) {
		if (a == 0) {
			return {0};
		}
		std::vector<T> res = {a / b};
		while (a != 0) {
			std::swap(a, b);
			res.emplace_back(a / b);
			a %= b;
		}
		return res;
	}
	template<i128::unsigned_integral T> std::vector<std::pair<T, T>> convergent(T a, T b) {
		std::vector<T> cf = cont_frac(a, b);
		std::vector<std::pair<T, T>> ret(cf.size());
		T p1 = 1, q1 = 0, p2 = 0, q2 = 1;
		for (u32 i = 0; i < cf.size(); ++i) {
			T p0 = cf[i] * p1 + p2, q0 = cf[i] * q1 + q2;
			ret[i] = std::make_pair(p0, q0);
			p2 = p1, q2 = q1;
			p1 = p0, q1 = q0;
		}
		return ret;
	}
}