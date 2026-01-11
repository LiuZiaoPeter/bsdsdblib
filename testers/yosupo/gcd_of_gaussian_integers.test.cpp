#define PROBLEM "https://judge.yosupo.jp/problem/gcd_of_gaussian_integers"

#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/gauss_int.hpp"

using GI = numtheo::GaussInt<i64>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		i32 a, b, c, d;
		std::cin >> a >> b >> c >> d;
		(numtheo::gauss_gcd(GI(a, b), GI(c, d))).outp_onlyspace();
		std::cout << '\n';
	}
	return 0;
}