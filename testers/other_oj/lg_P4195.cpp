#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/modint.hpp"
#include "../../numtheo/dis_log.hpp"

using MI = numtheo::ModInt32<-1>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 _a, p, _b;
	MI a, b;
	while (true) {
		std::cin >> _a >> p >> _b;
		if (!p) {
			break;
		}
		MI::set_mod(p);
		a = _a, b = _b;
		if (p == 1 || b == 1) {
			std::cout << "0\n";
			continue;
		}
		auto ans = numtheo::dis_log(a, b);
		if (ans.has_value()) {
			std::cout << ans.value() << '\n';
		} else {
			std::cout << "No Solution\n";
		}
	}
	return 0;
}