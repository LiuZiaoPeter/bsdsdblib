#include "../../numtheo/modint.hpp"
#include "../../numtheo/mul_inv.hpp"

using mip = numtheo::ModIntPr32<998244353>;

void init(int) {
	numtheo::O1inv_preproc<998244353, false>();
}

int inv(int x) {
	return static_cast<int>(numtheo::inv(mip(x, false)).value());
}