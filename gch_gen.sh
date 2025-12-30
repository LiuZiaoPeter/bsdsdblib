#!/bin/bash

shopt -s expand_aliases
alias compile_header='g++ -Ofast -std=c++20 -x c++-header 2>/dev/null'

if [ $# -eq 0 ] || [ "$1" == "--all-only" ]; then
	compile_header all.hpp
fi
if [ $# -eq 0 ] || [ "$1" == "--no-all" ]; then
	compile_header basics.hpp
	compile_header numtheo/euler_sieve.hpp
	compile_header numtheo/excrt.hpp
	compile_header numtheo/modint.hpp
	compile_header numtheo/o1gcd.hpp
	compile_header numtheo/pollard_rho.hpp
	compile_header numtheo/prim_root.hpp
	compile_header numtheo/prod_funcs.hpp
fi