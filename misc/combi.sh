#!/bin/bash
for arg in "$@"
do base=${arg##*/}
cat $base.fasta $base.err.fasta > $base.fasta.combined
done
