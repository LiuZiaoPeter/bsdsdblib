#!/bin/bash

./gch_gen.sh
exit_code=$?
find . -name "*.gch" -exec rm -f {} \;
exit $exit_code
