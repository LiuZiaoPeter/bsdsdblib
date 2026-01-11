#define PROBLEM "https://judge.yosupo.jp/problem/two_square_sum"

#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/square_decomp_2.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 q;
	std::cin >> q;
	while (q--) {
		u64 n;
		std::cin >> n;
		auto ans = numtheo::sqdecomp2_all(n);
		std::cout << ans.size() << '\n';
		for (auto i : ans) {
			std::cout << i.first << ' ' << i.second << '\n';
		}
	}
	return 0;
}
