#include <algorithm>
#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/pollard_rho.hpp"

u64 ans(u64 x) {
	if (numtheo::miller_rabin(x)) {
		return x;
	}
	u64 d = numtheo::pollard_rho(x);
	return std::max(ans(d), ans(x / d));
}

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u64 n;
		std::cin >> n;
		if (numtheo::miller_rabin(n)) {
			std::cout << "Prime\n";
		} else {
			std::cout << ans(n) << '\n';
		}
	}
	return 0;
}
