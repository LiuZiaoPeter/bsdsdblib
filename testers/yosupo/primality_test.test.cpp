#define PROBLEM "https://judge.yosupo.jp/problem/primality_test"

#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/pollard_rho.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u64 n;
		std::cin >> n;
		if (numtheo::miller_rabin(n)) {
			std::cout << "Yes\n";
		} else {
			std::cout << "No\n";
		}
	}
	return 0;
}