#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/prim_root.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u32 n, d;
		std::cin >> n >> d;
		auto ans = numtheo::prim_roots(n);
		std::cout << ans.size() << '\n';
		for (u32 i = d - 1; i < ans.size(); i += d) {
			std::cout << ans[i] << ' ';
		}
		std::cout << '\n';
	}
	return 0;
}
