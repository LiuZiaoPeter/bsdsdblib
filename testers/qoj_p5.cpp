#include "../numtheo/modint.hpp"
#include "../numtheo/mul_inv.hpp"

using mip = numtheo::ModIntPr<998244353, true>;

void init(int) {
	numtheo::O1inv_preproc<998244353, true>();
}

int inv(int x) {
	return static_cast<int>(numtheo::inv(mip(x, false)).value());
}