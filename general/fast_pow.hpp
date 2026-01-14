#pragma once

#include <bit>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "../basics.hpp"

template<class T, i128::unsigned_integral U> T qpow(T x, U y, const T &mul_iden = 1) {
	T ret = mul_iden;
	while (y) {
		if (y & 1) {
			ret *= x;
		}
		x *= x;
		y >>= 1;
	}
	return ret;
}

template<class T> class O1pow {
private:
	u32 B_log;
	std::vector<T> table_bs, table_gs;
public:
	O1pow() : B_log(0), table_bs(), table_gs() {}
	O1pow(T x, u64 lim, T mul_iden = 1) {
		B_log = std::bit_width(static_cast<u32>(std::sqrt(lim)) + 1);
		table_bs.resize(1u << B_log);
		T x_pow = mul_iden;
		for (u32 i = 0; !(i >> B_log); ++i, x_pow *= x) {
			table_bs[i] = x_pow;
		}
		table_gs.resize(1u << B_log);
		T x_to_B = x_pow, xB_pow = mul_iden;
		for (u32 i = 0; !(i >> B_log); ++i, xB_pow *= x_to_B) {
			table_gs[i] = xB_pow;
		}
	}
	template<i128::unsigned_integral U> T operator()(U x) {
		return table_bs[x & ((1u << B_log) - 1)] * table_gs[x >> B_log];
	}
};
