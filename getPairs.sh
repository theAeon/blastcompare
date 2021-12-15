#!/bin/bash
#export NCBI_API_KEY='<API KEY>'
#export EMAIL=<email>
find ./ -maxdepth 1 -type f -exec awk '!seen[$1]++{if ($3 >= 85) print > "{}best"}' {} \;
find ./*best -maxdepth 1 -type f -exec awk '{print $1 > "{}names"}' {} \;
find ./*best -maxdepth 1 -type f -exec awk '{print $2 > "{}IDs"}' {} \;
for file in ./*bestIDs; do
while read -r -u 9 line; do
blastdbcmd -db ../myblast/combined -entry "$line" >> "$file".fsa || efetch -db protein -id "$line" -format fasta >> "$file".fsa
done 9<"$file"
done

find ./chicken*names -maxdepth 1 -type f -exec seqkit grep -f {} -o {}.fsa ../input/chickendedup.fsa \;
find ./fish*names -maxdepth 1 -type f -exec seqkit grep -f {} -o {}.fsa ../input/fishdedup.fsa \;
find ./frog*names -maxdepth 1 -type f -exec seqkit grep -f {} -o {}.fsa ../input/frogdedup.fsa \;
find ./human*names -maxdepth 1 -type f -exec seqkit grep -f {} -o {}.fsa ../input/humandedup.fsa \;