#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/square_decomp_2.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u64 n;
		std::cin >> n;
		auto res = numtheo::sqdecomp2_all(n);
		std::cout << res.size() << '\n';
		for (auto [a, b] : res) {
			std::cout << a << ' ' << b << '\n';
		}
	}
	return 0;
}
