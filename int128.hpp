#pragma once

#include <algorithm>
#include <cstring>
#include <iostream>
#include <type_traits>

#include "basics.hpp"

// io
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
}