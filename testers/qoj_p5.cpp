#include "../numtheo/modint.hpp"
#include "../numtheo/mul_inv.hpp"

using mip = numtheo_n::ModIntPr<998244353, true>;

void init(int) {
	numtheo_n::O1inv_preproc<998244353, true>();
}

int inv(int x) {
	return static_cast<int>(numtheo_n::inv(mip(x, false)).value());
}