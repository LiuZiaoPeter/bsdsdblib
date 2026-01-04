#!/bin/bash

shopt -s expand_aliases
alias compile_cmd='g++ -Ofast -std=c++20 -x c++-header 2>/dev/null'

compile() {
	echo -n compiling $1 ""
	compile_cmd $1
	if [ $? -ne 0 ]; then 
		echo error
		exit
	else
		echo done
	fi
}

if [ $# -eq 0 ] || [ "$1" == "--all-only" ]; then
	compile all.hpp
fi
if [ $# -eq 0 ] || [ "$1" == "--no-all" ]; then
	while read line || [ -n "$line" ]; do
		if [ ${#line} -le 9 ] || [ "${line: 9: 1}" != "\"" ]; then
			continue
		fi
		cur_h="${line#*\"}"
		cur_h="${cur_h%\"*}"
		compile $cur_h
	done < all.hpp
fi