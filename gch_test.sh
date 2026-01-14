#!/bin/bash

./gch_gen.sh $1
exit_code=$?
find . -name "*.gch" -exec rm -f {} \;
exit $exit_code
