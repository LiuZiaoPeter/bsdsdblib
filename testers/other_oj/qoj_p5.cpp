#ifdef ONLINE_JUDGE
#include "inv.h"
#endif

#include "../../numtheo/modint.hpp"
#include "../../numtheo/farey_tech.hpp"

using mip = numtheo::ModIntPr32<998244353>;

void init(int) {
	numtheo::O1inv_preproc<998244353, false>();
}

int inv(int x) {
	return static_cast<int>(numtheo::inv(mip(x, false)).value());
}

#ifndef ONLINE_JUDGE

#include <iostream>
int main() {
	init(998244353);
	std::cout << inv(10000000) << std::endl; // 61689804
	return 0;
}

#endif
