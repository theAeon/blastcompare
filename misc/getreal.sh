#!/bin/bash
#export NCBI_API_KEY='<API KEY>'
#export EMAIL=<email>
lastresult=''
lastlookup=''
newColumn=''

while getopts ":s:f:" opt; do
  case $opt in
    s)
      species=$OPTARG
      echo "Querying Entrez for species $OPTARG" >&2
      ;;
    f)
      infile=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument" >&2
      exit 1
      ;;
  esac
done
testvar="$(awk '{gsub(/_/, " ", $1); print $1}' "$infile")"
totalLine="$(echo "$testvar" | wc -l)" 
counter=1
while read -r -u 9 line
do if [ "$lastlookup" = "$line" ]
      then
        newColumn+="$lastresult\n"
      else
        sanitized=$(echo "$line" | sed 'y/{[]}/{()}/')
        lastresult=$(esearch -db protein -query "$sanitized AND $species [ORGN]" | efetch -format acc | head -n 1)
        if [ "$lastresult" = "" ]
        then lastresult=sanitized
        fi
        newColumn+="$lastresult\n"
        lastlookup="$line"
     fi
   counterpercent=$(echo "scale=2; $counter / $totalLine * 100" | bc)
   echo -ne "${counterpercent%%.*}%\r"
   counter=$((counter +1))
done 9< <(echo "$testvar")
awk 'NR==FNR{a[NR]=$0;next} {$1=a[FNR]}1' <(echo -e "$newColumn") "$infile" > "$infile"2
