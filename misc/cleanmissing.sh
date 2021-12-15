#!/bin/bash
file=$1
sed 'y/{ }/{\t}/' "$file"2 | awk -F'\t' -v file="$file" '{if ($1 == "sanitized") {getline < file; print}}' > "$file"missing
temp=$(sed 'y/{ }/{\t}/' "$file"2 | awk -F'\t' '{if ($1 != "sanitized"){print}}')
rm "$file"2
echo "$temp" > "$file"2