#!/bin/bash
for arg in "$@"; do
while read -u 9 line; do
  seqkit rmdup -i $line -o $line.dedup
done 9<$arg
done
