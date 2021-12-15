#!/bin/bash
for arg in "$@"
do seqkit faidx --ignore-case -j 8 -o $arg.fai.gz $arg
done
