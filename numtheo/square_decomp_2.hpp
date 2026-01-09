// tested by lg_P2063

#pragma once

#include <complex>
#include <optional>
#include <utility>
#include <vector>

#include "../basics.hpp"
#include "../numtheo/cont_frac.hpp"
#include "../numtheo/pollard_rho.hpp"
#include "../numtheo/quad_res.hpp"

namespace numtheo {
	template<i128::liftable_unsigned T> std::pair<T, T> sqdecomp2_m4e1(T p) {
		using suT = i128::up_t<i128::make_signed_t<T>>;
		using MIP = numtheo::ModIntPr<modint_inner + 6, std::is_same_v<T, u64>>; // mod : p
		MIP::set_mod(p);
		T c = numtheo::sqrt(MIP(p - 1, false)).value().value();
		auto conv = convergent(c, p);
		auto it = std::prev(conv.end());
		while (static_cast<i128::up_t<T>>(it->second) * it->second > p) {
			--it;
		}
		auto [a, b] = *it;
		return {std::abs(static_cast<i128::make_signed_t<T>>(static_cast<suT>(b * c) - static_cast<suT>(p * a))), b};
	}
	template<i128::liftable_unsigned T> std::optional<std::pair<T, T>> sqdecomp2_one(T x) {
		using cplx = std::complex<i128::make_signed_t<T>>;
		auto pf = numtheo::prime_factors(x);
		cplx ret(1, 0);
		for (auto [p, a] : pf) {
			if (p == 2) {
				ret *= qpow(cplx(1, 1), a);
			} else if ((p & 3) == 3) {
				if (a & 1) {
					return std::nullopt;
				}
				ret *= qpow(p, a >> 1);
			} else {
				std::pair<T, T> dec = sqdecomp2_m4e1(p);
				cplx dec_c(dec.first, dec.second);
				ret *= qpow(dec_c, a);
			}
		}
		return std::make_pair(std::abs(ret.real()), std::abs(ret.imag()));
	}
	template<i128::liftable_unsigned T> std::vector<std::pair<T, T>> sqdecomp2_all(T x) {
		using cplx = std::complex<i128::make_signed_t<T>>;
		auto pf = numtheo::prime_factors(x);
		std::vector<cplx> r1 = {cplx(1, 0)};
		for (auto [p, a] : pf) {
			if (p == 2) {
				cplx tmp = qpow(cplx(1, 1), a);
				for (auto &z : r1) {
					z *= tmp;
				}
			} else if ((p & 3) == 3) {
				if (a & 1) {
					return {};
				}
				T tmp = qpow(p, a >> 1);
				for (auto &z : r1) {
					z *= tmp;
				}
			} else {
				std::pair<T, T> dec = sqdecomp2_m4e1(p);
				cplx dec_c(dec.first, dec.second);
				std::vector<cplx> pws = {cplx(1, 0)};
				pws.resize(a + 1);
				for (u32 i = 1; i <= a; ++i) {
					pws[i] = pws[i - 1] * dec_c;
				}
				std::vector<cplx> r2;
				r2.reserve(r1.size() * (a + 1));
				for (auto z : r1) {
					for (u32 i = 0; i <= a; ++i) {
						r2.emplace_back(z * pws[i] * std::conj(pws[a - i]));
					}
				}
				r1 = r2;
			}
		}
		std::vector<std::pair<T, T>> ret;
		for (auto z : r1) {
			ret.emplace_back(std::abs(z.real()), std::abs(z.imag()));
			ret.emplace_back(std::abs(z.imag()), std::abs(z.real()));
		}
		std::sort(ret.begin(), ret.end());
		ret.erase(std::unique(ret.begin(), ret.end()), ret.end());
		return ret;
	}
}