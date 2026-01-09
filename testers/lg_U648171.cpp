#include <iostream>

#include "../basics.hpp"
#include "../numtheo/modint.hpp"
#include "../numtheo/dis_log.hpp"

using MI = numtheo::ModInt64<-1>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	bool t;
	u64 p;
	u32 n;
	std::cin >> t >> p >> n;
	MI::set_mod(p);
	MI a;
	std::cin >> a;
	std::vector<MI> b(n);
	for (MI &i : b) {
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