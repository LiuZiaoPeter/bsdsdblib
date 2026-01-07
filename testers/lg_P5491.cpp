#include <iostream>
#include <vector>

#include "../basics.hpp"
#include "../numtheo/modint.hpp"
#include "../numtheo/quad_res.hpp"

using mip = numtheo_n::ModIntPr<-1, false>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u32 p, n;
		std::cin >> n >> p;
		if (n == 0) {
			std::cout << "0\n";
			continue;
		}
		mip::set_mod(p);
		auto ans = numtheo_n::sqrt(mip(n));
		if (ans.has_value() == false) {
			std::cout << "Hola!\n";
		} else {
			std::cout << ans.value().value() << ' ' << p - ans.value().value() << '\n';
		}
	}
	return 0;
}