#pragma once

#include "../basics.hpp"
#include "../int128.hpp"

namespace numtheo_n {

	template<class Derived, i64 P, bool _64 = false> class ModIntBase {
		/*
		lt 0 for dynamic, le -1073741824 for internal use
		occupied P shown below:
		-1073741824 : template<i32 P> ModInt<P> std::optional<u32> ord(ModInt<P>)
		*/
	public:
		using val_t = std::conditional_t<_64, u64, u32>;
		using mul_t = std::conditional_t<_64, __uint128_t, u64>;
	protected:
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
		friend std::istream &operator>>(std::istream &in, Derived &x) {
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
	
	template<i64 P, bool _64> class ModIntPr : public ModIntBase<ModIntPr<P, _64>, P, _64> { // P prime
	private:
		using Base = ModIntBase<ModIntPr<P, _64>, P, _64>;
		using Base::Base;
	public:
		using Base::val_t;
		using Base::mul_t;
		using Base::mod;
		using Base::set_mod;
	};

	template<i64 P, bool _64> class ModInt : public ModIntBase<ModInt<P, _64>, P, _64> {
	private:
		using Base = ModIntBase<ModInt<P, _64>, P, _64>;
		using Base::Base;
	public:
		using Base::val_t;
		using Base::mul_t;
		using Base::mod;
		using Base::set_mod;
	};

}