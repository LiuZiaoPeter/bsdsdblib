#define PROBLEM "https://judge.yosupo.jp/problem/associative_array"

#include <iostream>

#include "../../ds/hash_table.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 q;
	std::cin >> q;
	static ds::hash_table<u64, u64> t;
	while (q--) {
		u32 o;
		u64 k, v;
		std::cin >> o;
		if (o == 0) {
			std::cin >> k >> v;
			t[k] = v;
		} else {
			std::cin >> k;
			std::cout << t[k] << '\n';
		}
	}
	return 0;
}
