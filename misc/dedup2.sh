#!/bin/bash
export NCBI_API_KEY='c9838a563f0e8d1240d2eb6d62afee466908'
for arg in "$@"; do
while read -u 9 line; do
  seqkit rmdup -i $line -o $line.dedup
done 9<$arg
done
