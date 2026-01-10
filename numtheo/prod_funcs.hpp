#pragma once

#include "../basics.hpp"
#include "pollard_rho.hpp"

namespace numtheo {
	template<i128::liftable_unsigned T> T phi(T x) {
		auto prf = prime_factors(x);
		for (auto i : prf) {
			x = x / i.first * (i.first - 1);
		}
		return x;
	}
}
