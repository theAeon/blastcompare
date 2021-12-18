#!/bin/bash
for arg in "$@"; do
while read -u 9 line; do
  ID=$(awk -F'Skipped ' '{print $2}' <<< $line)
  esearch -db protein -query $ID |\
    efetch -format fasta >> $arg.fasta
done 9<$arg
done
