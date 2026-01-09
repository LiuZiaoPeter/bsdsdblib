#include <iostream>

#include "../basics.hpp"
#include "../numtheo/modint.hpp"
#include "../numtheo/dis_log.hpp"

using mi = numtheo::ModInt<-1, true>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	bool t;
	u64 p;
	u32 n;
	std::cin >> t >> p >> n;
	mi::set_mod(p);
	mi a;
	std::cin >> a;
	std::vector<mi> b(n);
	for (mi &i : b) {
		std::cin >> i;
	}
	auto ans = numtheo::dis_logs(a, b);
	for (auto i : ans) {
		if (i.has_value()) {
			std::cout << i.value() << ' ';
		} else {
			std::cout << "-1 ";
		}
	}
	std::cout << '\n';
	return 0;
}