#include <iostream>
#include <vector>

#include "../numtheo/modint.hpp"
#include "../numtheo/mul_inv.hpp"
#include "../numtheo/dis_log.hpp"
#include "../numtheo/fast_dis_ln.hpp"
#include "../numtheo/quad_res.hpp"

constexpr u64 modp = static_cast<u64>(1e10) + 19;
using mip = numtheo_n::ModIntPr<modp, true>;
const mip g = 2;
constexpr u64 mod = static_cast<u64>(1e10);
using mi = numtheo_n::ModInt<mod, true>;

int main() {
	// inv 6568999657
	std::cout << numtheo_n::inv(mip(998244353)).value() << std::endl;
	numtheo_n::O1inv_preproc<modp, true>();
	std::cout << numtheo_n::inv(mip(998244353)).value() << std::endl;
	std::cout << (mip(1) / mip(998244353)).value() << std::endl;
	// dis log
	std::vector<mip> qp = {mip(2), mip(4), mip(8), mip(9286000679)};
	for (mip i : qp) {
		std::cout << numtheo_n::dis_log(g, i).value() << ' ';
	}
	std::cout << std::endl;
	for (auto i : numtheo_n::dis_logs(g, qp)) {
		std::cout << i.value() << ' ';
	}
	std::cout << std::endl;
	numtheo_n::dis_ln_preproc(g);
	for (mip i : qp) {
		std::cout << numtheo_n::fast_dis_ln(i) << ' ';
	}
	std::cout << std::endl;
	std::vector<mi> q = {mi(2), mi(4), mi(8), mi(6703205376)};
	for (mi i : q) {
		std::cout << numtheo_n::dis_log(mi(2), i).value() << ' ';
	}
	std::cout << std::endl;
	for (auto i : numtheo_n::dis_logs(mi(2), q)) {
		std::cout << i.value() << ' ';
	}
	std::cout << std::endl;
	// quad res 998244353
	std::cout << numtheo_n::sqrt(mip(6403054227)).value().value() << std::endl;
	return 0;
}