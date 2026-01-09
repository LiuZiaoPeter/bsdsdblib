#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <iostream>
#include <type_traits>

using u8 = uint8_t;
using i8 = int8_t;
using u16 = uint16_t;
using i16 = int16_t;
using u32 = uint32_t;
using i32 = int32_t;
using u64 = uint64_t;
using i64 = int64_t;
using __u128 = __uint128_t;
using __i128 = __int128_t;

std::ostream &operator<<(std::ostream &out, __uint128_t x) {
	if (x == 0) {
		out << "0";
		return out;
	}
	std::string _o;
	while (x != 0) {
		_o += char((x % 10) | 48);
		x /= 10;
	}
	std::reverse(_o.begin(), _o.end());
	out << _o;
	return out;
}
std::ostream &operator<<(std::ostream &out, __int128_t x) {
	if (x < 0) {
		out << "-";
		x *= -1;
	}
	out << static_cast<__uint128_t>(x);
	return out;
}

namespace i128 {	
	// concepts
	template<class T> concept unsigned_integral = std::unsigned_integral<T> || std::is_same_v<T, __uint128_t>;
	template<class T> concept signed_integral = std::signed_integral<T> || std::is_same_v<T, __int128_t>;
	template<class T> concept liftable_unsigned = std::unsigned_integral<T>;
	template<class T> concept liftable_signed = std::signed_integral<T>;
	// type traits
	template<class T> struct make_unsigned {
		using type = std::make_unsigned_t<T>;
	};
	template<> struct make_unsigned<__int128_t> {
		using type = __uint128_t;
	};
	template<class T> using make_unsigned_t = make_unsigned<T>::type;
	template<class T> struct make_signed {
		using type = std::make_signed_t<T>;
	};
	template<> struct make_signed<__uint128_t> {
		using type = __int128_t;
	};
	template<class T> using make_signed_t = make_signed<T>::type;
	template<class T> struct up {
		using type = void;
	};
	template<> struct up<u8> {
		using type = u16;
	};
	template<> struct up<u16> {
		using type = u32;
	};
	template<> struct up<u32> {
		using type = u64;
	};
	template<> struct up<u64> {
		using type = __u128;
	};
	template<> struct up<i8> {
		using type = i16;
	};
	template<> struct up<i16> {
		using type = i32;
	};
	template<> struct up<i32> {
		using type = i64;
	};
	template<> struct up<i64> {
		using type = __i128;
	};
	template<class T> using up_t = up<T>::type;
}

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

#define __func_str__ static_cast<std::string>(__func__)