#define PROBLEM "https://judge.yosupo.jp/problem/discrete_logarithm_mod"

#include <iostream>
#include <optional>

#include "../../basics.hpp"
#include "../../numtheo/modint.hpp"
#include "../../numtheo/dis_log.hpp"

using MI = numtheo::ModInt32<-1>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u32 x, y, m;
		std::cin >> x >> y >> m;
		if (m == 1) {
			std::cout << "0\n";
			continue;
		}
		if (x == 0) {
			if (y == 0) {
				std::cout << "1\n";
			} else if (y == 1) {
				std::cout << "0\n";
			} else {
				std::cout << "-1\n";
			}
			continue;
		}
		if (y == 1) {
			std::cout << "0\n";
			continue;
		}
		MI::set_mod(m);
		auto ans = numtheo::dis_log(MI(x, false), MI(y, false));
		if (ans.has_value() == true) {
			std::cout << ans.value() << '\n';
		} else {
			std::cout << "-1\n";
		}
	}
	return 0;
}
