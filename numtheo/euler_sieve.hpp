#pragma once

#include <numeric>
#include <vector>

#include "../basics.hpp"

namespace numtheo {
	std::vector<u32> primes, mpf;
	u32 enumerated_prime;
	void euler_sieve(u32 N) {
		if (mpf.size() > N) {
			return;
		}
		primes.clear();
		mpf.resize(N + 1);
		std::iota(mpf.begin(), mpf.end(), 0);
		enumerated_prime = N;
		for (u32 i = 2; i <= N; ++i) {
			if (mpf[i] == i) {
				primes.emplace_back(i);
			}
			for (u32 j : primes) {
				u64 k = static_cast<u64>(i) * j;
				if (k > N) {
					break;
				}
				mpf[k] = j;
				if (i % j == 0) {
					break;
				}
			}
		}
	}
	void enum_prime(u32 N) {
		if (enumerated_prime >= N) {
			return;
		}
		primes.clear();
		std::vector<bool> is_prime(N + 1, 1);
		enumerated_prime = N;
		for (u32 i = 2; i <= N; ++i) {
			if (is_prime[i] == true) {
				primes.emplace_back(i);
			}
			for (u32 j : primes) {
				u64 k = static_cast<u64>(i) * j;
				if (k > N) {
					break;
				}
				is_prime[k] = false;
				if (i % j == 0) {
					break;
				}
			}
		}
	}
}
