#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_mod"

#include <iostream>
#include <optional>

#include "../../basics.hpp"
#include "../../numtheo/modint.hpp"
#include "../../numtheo/quad_res.hpp"

using MIP = numtheo::ModIntPr32<-1>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u32 y, p;
		std::cin >> y >> p;
		MIP::set_mod(p);
		auto ans = numtheo::sqrt(MIP(y, false));
		if (ans.has_value() == true) {
			std::cout << ans.value().value() << '\n';
		} else {
			std::cout << "-1\n";
		}
	}
	return 0;
}
