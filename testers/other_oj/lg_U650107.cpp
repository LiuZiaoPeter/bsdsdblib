#include <iostream>

#include "../../basics.hpp"
#include "../../general/fast_pow.hpp"

int main() {
	u64 a;
	std::cin >> a;
	O1pow<u64> pw(a, (1ull << 46) - 1);
	u32 n;
	std::cin >> n;
	u64 b_last, seed, ans = 0;
	std::cin >> b_last >> seed;
	for (u32 i = 1; i <= n; ++i) {
		u64 b_cur = (pw(b_last) ^ (seed + i)) & ((1ull << 46) - 1);
		ans += i * pw(b_cur);
		b_last = b_cur;
	}
	std::cout << ans << std::endl;
	return 0;
}
