#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/modint.hpp"
#include "../../numtheo/o1gcd.hpp"

using MIP = numtheo::ModIntPr32<998244353>;

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	numtheo::O1gcd_preproc(1000000);
	u32 n;
	std::cin >> n;
	std::vector<u32> a(n), b(n);
	for (u32 &i : a) {
		std::cin >> i;
	}
	for (u32 &i : b) {
		std::cin >> i;
	}
	for (u32 i = 0; i < n; ++i) {
		MIP ans = 0, coe = 1;
		for (u32 j = 0; j < n; ++j) {
			coe *= i + 1;
			ans += MIP(numtheo::O1gcd(a[i], b[j]), false) * coe;
		}
		std::cout << ans.value() << '\n';
	}
	return 0;
}