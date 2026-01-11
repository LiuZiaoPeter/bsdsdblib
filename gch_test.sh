#!/bin/bash

./gch_gen.sh
find . -name "*.gch" -exec rm -f {} \;