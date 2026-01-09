#include <iostream>

#include "../basics.hpp"
#include "../numtheo/modint.hpp"
#include "../numtheo/dis_log.hpp"

using mip = numtheo::ModIntPr<-1, false>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 p;
	std::cin >> p;
	mip::set_mod(p);
	mip b, n;
	std::cin >> b >> n;
	if (n == 1) {
		std::cout << "0\n";
	} else {
		auto ans = numtheo::dis_log(b, n);
		if (ans.has_value()) {
			std::cout << ans.value() << std::endl;
		} else {
			std::cout << "no solution\n";
		}
	}
	return 0;
}