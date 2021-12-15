""" Takes fasta sequences from two ordered fasta files, in which each pair is a known alignment.
Developed against Python 3.9.7, BioPython 1.78, and NumPy 1.19.5
Tested with Muscle v3.8.1551.
"""
from argparse import ArgumentParser
from os import get_terminal_size
from collections import Counter
from io import StringIO
from subprocess import run
from typing import Literal, Tuple
from Bio import AlignIO, SeqIO
term_size = get_terminal_size()
AminoAcid = Literal["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "-"]
Clustal = Literal['','.',':','*']
Position = Tuple[int, Clustal, AminoAcid, AminoAcid]
PairDiff = list[Position]

parser = ArgumentParser(prog="get_diff",
                        description="a script that takes paired fasta sequences from two species and returns count of mismatched AAs"
                        )
parser.add_argument("species", action="extend", nargs="+", help="species base used as input for getPairs.sh")
parsed = parser.parse_args()


class SpeciesCompare():
    def __init__(self,base: str) -> None:
        self.seq_ids=SeqIO.parse("%sbestIDs.fsa"  % (base), "fasta")
        self.seqnames=SeqIO.parse("%sbestnames.fsa"  % (base), "fasta")
        self.seqidlist=list(self.seq_ids)
        self.seqnamelist=list(self.seqnames)
    def getpairlist(self, i: int):
        combistr = format(self.seqnamelist[i], "fasta") + format(self.seqidlist[i], "fasta")
        clustal = run(
            ["muscle", "-quiet", "-clwstrict", "-maxiters", "1", "-diags", "-sv", "-distance1", "kbit20_3"],
            capture_output=True,
            text=True,
            input=combistr,
            check=True
            )
        alignment = AlignIO.read(StringIO(clustal.stdout), "clustal")
        pairlist: PairDiff = []
        for kndex, letter in enumerate(alignment.column_annotations['clustal_consensus']):
            base1 = alignment[0, kndex]
            base2 = alignment[1, kndex]
            diff = letter
            loc = kndex
            if diff != "*":
                pairlist.append((loc, diff, base1, base2))
        return pairlist

listComparisons: list[SpeciesCompare] = []
for comparator in parsed.species:
    local = SpeciesCompare(comparator)
    testAll=[local.getpairlist(i) for i in range(len(local.seqidlist))]
    flat=[item for sublist in testAll for item in sublist]
    testA=Counter(list(list(zip(*flat))[2]))
    testB=Counter(list(list(zip(*flat))[3]))
    print(comparator.split('v')[0])
    print(testA)
    print(comparator.split('v')[1])
    print(testB)
    print(u'\u2500' * term_size.columns)
