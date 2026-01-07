#include <iostream>

#include "../basics.hpp"
#include "../numtheo/euler_sieve.hpp"

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr), std::cout.tie(nullptr);
	u32 n, q;
	std::cin >> n >> q;
	numtheo_n::euler_sieve(n);
	while (q--) {
		u32 x;
		std::cin >> x;
		std::cout << numtheo_n::primes[x - 1] << '\n';
	}
	return 0;
}