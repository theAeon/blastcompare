#!/bin/bash
for arg in "$@"
do base=${arg##*/}
mv $base.fasta.combined.dedup $base.fasta
done
