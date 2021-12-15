# blastcompare

This repository contains a series of work-in-progress scripts for the purpose of comparing amino acid frequency across species.
The misc folder contains files not needed for the primary pipeline.

## Usage
Run a blast search as in runBlast.sh with the full proteome of one species against another. Ensure output contains identity % in the third colum for filtering these files should be named <query>v<target>. Use getPairs.sh with e-utilities and blastdbcmd to get paired fasta files in the output directory, and use get_diff.py to print amino acid frequencies. More interactivity coming soon-maybe.

## Dependencies (ish)
Python 3.9, Biopython 1.78, MUSCLE 3.8, a recent copy of BLAST+ and a recent copy of NCBI EUtils.
