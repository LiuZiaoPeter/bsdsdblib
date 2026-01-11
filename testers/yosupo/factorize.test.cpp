#define PROBLEM "https://judge.yosupo.jp/problem/factorize"

#include <iostream>
#include <vector>

#include "../../basics.hpp"
#include "../../numtheo/pollard_rho.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u64 n;
		std::cin >> n;
		std::vector<u64> ans;
		numtheo::prime_factors(n, ans);
		std::sort(ans.begin(), ans.end());
		std::cout << ans.size() << ' ';
		for (u64 i : ans) {
			std::cout << i << ' ';
		}
		std::cout << '\n';
	}
	return 0;
}
