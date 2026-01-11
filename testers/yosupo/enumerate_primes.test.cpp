#define PROBLEM "https://judge.yosupo.jp/problem/enumerate_primes"

#include <iostream>

#include "../../basics.hpp"
#include "../../numtheo/euler_sieve.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 n, a, b;
	std::cin >> n >> a >> b;
	numtheo::enum_prime(n);
	std::cout << numtheo::primes.size() << ' ' << (numtheo::primes.size() - 1 - b) / a + 1 << std::endl;
	for (u32 i = b; i < numtheo::primes.size(); i += a) {
		std::cout << numtheo::primes[i] << ' ';
	}
	std::cout << std::endl;
	return 0;
}
