#include <iostream>

#include "../basics.hpp"
#include "../numtheo/modint.hpp"
#include "../numtheo/dis_log.hpp"

using mi = numtheo::ModInt<-1, true>;

int main() {
	u64 p;
	std::cin >> p;
	mi::set_mod(p);
	mi a;
	std::cin >> a;
	auto ans = numtheo::ord(a);
	if (ans.has_value()) {
		std::cout << ans.value() << std::endl;
	} else {
		std::cout << "-1\n";
	}
	return 0;
}