#define PROBLEM "https://judge.yosupo.jp/problem/primitive_root"

#include <iostream>
#include <optional>

#include "../../basics.hpp"
#include "../../numtheo/prim_root.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 T;
	std::cin >> T;
	while (T--) {
		u64 n;
		std::cin >> n;
		std::cout << numtheo::min_prim_root(n).value() << '\n';
	}
	return 0;
}
