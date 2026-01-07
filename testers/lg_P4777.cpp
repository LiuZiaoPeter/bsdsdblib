#include <iostream>
#include <vector>

#include "../basics.hpp"
#include "../numtheo/excrt.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 n;
	std::cin >> n;
	std::vector<u64> p(n), a(n);
	for (u32 i = 0; i < n; ++i) {
		std::cin >> p[i] >> a[i];
	}
	auto ans = numtheo_n::excrt(a, p);
	std::cout << ans.value().first << std::endl;
	return 0;
}