#include <iostream>
#include <vector>

#include "../basics.hpp"
#include "../numtheo/modint.hpp"
#include "../numtheo/fast_dis_ln.hpp"

using mip = numtheo::ModIntPr<-1, false>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 p;
	std::cin >> p;
	mip::set_mod(p);
	mip g;
	std::cin >> g;
	numtheo::dis_ln_preproc(g);
	u32 q;
	std::cin >> q;
	while (q--) {
		mip x;
		std::cin >> x;
		std::cout << numtheo::fast_dis_ln(x) << '\n';
	}
	return 0;
}