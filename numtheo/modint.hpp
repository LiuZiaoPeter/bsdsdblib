#pragma once

#include <numeric>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

#include "../basics.hpp"
#include "../int128.hpp"
#include "euler_sieve.hpp"
#include "pollard_rho.hpp"
#include "prod_funcs.hpp"

namespace numtheo_n {
	template<i64 P, bool _64 = false> class MIP;
	template<i32 P> class MI;
	// multiplicative inverse
	template<i128::signed_integral T, i64 P, bool _64 = false> MIP<P, _64> qpow_signed(MIP<P, _64>, T);
	template<i64 P, bool _64 = false> MIP<P, _64> inv(MIP<P, _64>);
	// discrete log
	template<i64 P, bool _64 = false> std::optional<std::conditional_t<_64, u64, u32>> dis_log(MIP<P, _64>, MIP<P, _64>);
	template<i64 P, bool _64 = false> std::vector<std::optional<std::conditional_t<_64, u64, u32>>> dis_logs(MIP<P, _64>, std::vector<MIP<P, _64>>);
	template<i64 P, bool _64 = false> std::conditional_t<_64, u64, u32> fast_dis_ln(MIP<P, _64>);
	template<i64 P, bool _64 = false> std::conditional_t<_64, u64, u32> ord(MIP<P, _64>);
	template<i32 P> std::optional<u32> dis_log(MI<P>, MI<P>);
	template<i32 P> std::vector<std::optional<u32>> dis_logs(MI<P>, std::vector<MI<P>>);
	template<i32 P> std::optional<u32> ord(MI<P> x);
	// quadradic residue
	template<i64 P, bool _64 = false> i32 legendre(MIP<P, _64>);
	template<i64 P, bool _64 = false> std::optional<MIP<P, _64>> sqrt(MIP<P, _64>);
	template<class Derived, i64 P, bool _64 = false> class ModIntBase { // lt 0 for dynamic, le -1073741824 for internal use
		/* 
		occupied P shown below:
		-1073741824 : template<i32 P> MI<P> std::optional<u32> ord(MI<P>)
		*/
	protected:
		using val_t = std::conditional_t<_64, u64, u32>;
		using mul_t = std::conditional_t<_64, __uint128_t, u64>;
		val_t val;
		inline static val_t dyn_mod = 0;
	public:
		static void set_mod(val_t p) {
			dyn_mod = p;
		}
		static val_t mod() {
			if (P < 0) {
				return dyn_mod;
			}
			return P;
		}
		ModIntBase() : val(0) {}
		template<i128::unsigned_integral T> ModIntBase(T v, bool need_mod = true) : val(
			static_cast<val_t>(need_mod ? v % mod() : v)
		) {}
		template<i128::signed_integral T> ModIntBase(T v, bool need_mod = true) : val(
			static_cast<val_t>(need_mod ? (v < 0 ? mod() - (-v) % mod() : v % mod()) : v)
		) {}
		void read_by_mod() {
			std::string s;
			std::cin >> s;
			val = 0;
			for (char c : s) {
				val = (static_cast<mul_t>(val) * 10 + (c ^ 48)) % mod();
			}
		}
		friend std::istream &operator>>(std::istream &in, ModIntBase<Derived, P> &x) {
			in >> x.val;
			return in;
		}
		val_t value() const {
			return val;
		}
		Derived operator+(Derived x) const {
			Derived ret;
			ret.val = val + x.val;
			if (ret.val >= mod()) {
				ret.val -= mod();
			}
			return ret;
		}
		Derived &operator+=(Derived x) {
			val += x.val;
			if (val >= mod()) {
				val -= mod();
			}
			return static_cast<Derived&>(*this);
		}
		Derived operator-() const {
			Derived ret;
			ret.val = mod() - val;
			return ret;
		}
		Derived operator-(Derived x) const {
			Derived ret;
			ret.val = val + mod() - x.val;
			if (ret.val >= mod()) {
				ret.val -= mod();
			}
			return ret;
		}
		Derived &operator-=(Derived x) {
			val += mod() - x.val;
			if (val >= mod()) {
				val -= mod();
			}
			return static_cast<Derived&>(*this);
		}
		Derived operator*(Derived x) const {
			return static_cast<Derived>(static_cast<mul_t>(val) * x.val);
		}
		Derived &operator*=(Derived x) {
			val = static_cast<val_t>(static_cast<mul_t>(val) * x.val % mod());
			return static_cast<Derived&>(*this);
		}
		friend bool operator==(Derived x, Derived y) {
			return x.val == y.val;
		}
		friend bool operator!=(Derived x, Derived y) {
			return x.val != y.val;
		}
	};
	template<i64 P, bool _64> class MIP : public ModIntBase<MIP<P, _64>, P, _64> { // P prime
	private:
		using Base = ModIntBase<MIP<P, _64>, P, _64>;
		using Base::Base;
		using typename Base::val_t;
		using typename Base::mul_t;
		// O(1) inv
		inline static std::vector<MIP<P, _64>> inv_v;
		inline static val_t cbrtP_log2, cbrtP, cbrtP2;
		inline static bool O1inv_mode = false;
		inline static std::vector<std::pair<val_t, val_t>> farey_v, farey_prec, farey_succ;
		// fast discrete ln
		inline static std::vector<val_t> lesqrt_ln;
	public:
		using Base::mod;
		using Base::set_mod;
		// multiplicative inverse
		template<i128::signed_integral T> friend MIP<P, _64> qpow_signed(MIP<P, _64>, T);
		template<i128::unsigned_integral T> static void lin_inv_preproc(T);
		static void O1inv_preproc();
		friend MIP<P, _64> inv<>(MIP<P, _64>);
		MIP<P, _64> operator/(MIP<P, _64> x) const {
			return *this * inv(x);
		}
		MIP<P, _64> &operator/=(MIP<P, _64> x) {
			return *this = *this * inv(x);
		}
		// discrete log
		friend std::optional<val_t> dis_log<>(MIP<P, _64>, MIP<P, _64>);
		friend std::vector<std::optional<val_t>> dis_logs<>(MIP<P, _64>, std::vector<MIP<P, _64>>);
		static void dis_ln_preproc(MIP<P, _64>);
		friend val_t fast_dis_ln<>(MIP<P, _64>);
		friend val_t ord<>(MIP<P, _64>);
		// quadradic residue
		friend i32 legendre<>(MIP<P, _64>);
		friend std::optional<MIP<P, _64>> sqrt<>(MIP<P, _64>);
	};
	template<i32 P> class MI : public ModIntBase<MI<P>, P> {
	private:
		using Base = ModIntBase<MI<P>, P>;
		using Base::Base;
	public:
		using Base::mod;
		using Base::set_mod;
		// discrete log
		friend std::optional<u32> dis_log<>(MI<P>, MI<P>);
		friend std::vector<std::optional<u32>> dis_logs<>(MI<P>, std::vector<MI<P>>);
		friend std::optional<u32> ord<>(MI<P>);
	};
}

#include "modint_class/fast_dis_ln.hpp"
#include "modint_class/dis_log.hpp"
#include "modint_class/mul_inv.hpp"
#include "modint_class/quad_res.hpp"