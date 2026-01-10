#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/modint.hpp"
#include "../../numtheo/dis_log.hpp"

using MI = numtheo::ModInt64<-1>;

int main() {
	u64 p;
	std::cin >> p;
	MI::set_mod(p);
	MI a;
	std::cin >> a;
	auto ans = numtheo::ord(a);
	if (ans.has_value()) {
		std::cout << ans.value() << std::endl;
	} else {
		std::cout << "-1\n";
	}
	return 0;
}
