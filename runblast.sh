#!/bin/bash

THREADS=1
OUTFMT="6"
while getopts ":t:o:" opt; do
  case $opt in
    t)
      THREADS=$OPTARG
      echo "Running blast with $OPTARG threads." >&2
      ;;
    o)
      OUTFMT=$OPTARG
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

echo -ne 'Task 01/12\r'
blastp -db gallus.fsa -query input/humandedup.fsa -out output/humanvchicken -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 02/12\r'
blastp -db laevis.fsa -query input/humandedup.fsa -out output/humanvfrog -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 03/12\r'
blastp -db rerio.fsa -query input/humandedup.fsa -out output/humanvfish -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 04/12\r'
blastp -db sapiens.fsa -query input/fishdedup.fsa -out output/fishvhuman -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 05/12\r'
blastp -db laevis.fsa -query input/chickendedup.fsa -out output/chickenvfrog -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 06/12\r'
blastp -db rerio.fsa -query input/chickendedup.fsa -out output/chickenvfish -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 07/12\r'
blastp -db sapiens.fsa -query input/chickendedup.fsa -out output/chickenvhuman -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 08/12\r'
blastp -db laevis.fsa -query input/fishdedup.fsa -out output/fishvfrog -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 09/12\r'
blastp -db gallus.fsa -query input/fishdedup.fsa -out output/fishvchicken -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 10/12\r'
blastp -db gallus.fsa -query input/frogdedup.fsa -out output/frogvchicken -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 11/12\r'
blastp -db sapiens.fsa -query input/frogdedup.fsa -out output/frogvhuman -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90

echo -ne 'Task 12/12\r'
blastp -db rerio.fsa -query input/frogdedup.fsa -out output/frogvfish -num_threads $THREADS -outfmt "$OUTFMT" -mt_mode 1 -qcov_hsp_perc 90
