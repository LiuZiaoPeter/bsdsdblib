#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "../basics.hpp"
#include "euler_sieve.hpp"

namespace numtheo {
	std::vector<std::array<u32, 3>> gcd_decomp;
	std::vector<std::vector<u32>> gcd_table;
	void O1gcd_preproc(u32 N) {
		u32 sqrtN = static_cast<u32>(std::sqrt(N)) + 1;
		euler_sieve(N);
		gcd_decomp.clear(), gcd_decomp.resize(N + 1);
		for (u32 i = 1; i <= N; ++i) {
			if (mpf[i] == i) {
				gcd_decomp[i] = {1, 1, i};
			} else {
				gcd_decomp[i] = gcd_decomp[i / mpf[i]];
				gcd_decomp[i][0] *= mpf[i];
				std::sort(gcd_decomp[i].begin(), gcd_decomp[i].end());
			}
		}
		gcd_table.assign(sqrtN + 1, std::vector<u32>(sqrtN + 1, 0));
		for (u32 i = 0; i <= sqrtN; ++i) {
			for (u32 j = 0; j <= sqrtN; ++j) {
				if (std::min(i, j) == 0) {
					gcd_table[i][j] = std::max(i, j);
				} else if (std::min(i, j) == 1) {
					gcd_table[i][j] = 1;
				} else {
					u32 dg = gcd_table[i][j / mpf[j]];
					if ((i / dg) % mpf[j] == 0) {
						gcd_table[i][j] = dg * mpf[j];
					} else {
						gcd_table[i][j] = dg;
					}
				}
			}
		}
	}
	u32 O1gcd(u32 x, u32 y) {
		assure(x < gcd_decomp.size(), "x={} out of gcd preproc range", x);
		assure(y < gcd_decomp.size(), "y={} out of gcd preproc range", y);
		u32 ret = 1;
		for (u32 i : gcd_decomp[x]) {
			u32 g;
			if (i >= gcd_table.size()) {
				g = y % i ? 1 : i;
			} else {
				g = gcd_table[i][y % i];
			}
			ret *= g;
			y /= g;
		}
		return ret;
	}
}
