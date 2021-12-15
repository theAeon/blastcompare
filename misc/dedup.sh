#!/bin/bash
for arg in "$@"
do seqkit rmdup --ignore-case --by-seq -j 13 -o $arg.gz $arg
done
