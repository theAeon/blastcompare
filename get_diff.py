""" Takes fasta sequences from two ordered fasta files, in which each pair is a known alignment.
Developed against Python 3.9.7, BioPython 1.78, and NumPy 1.19.5
Tested with Muscle v3.8.1551.
"""
from argparse import ArgumentParser
from os import get_terminal_size
from collections import Counter
from pprint import pprint
from io import StringIO
from subprocess import run
#from typing import Literal, Tuple
#from types import new_class
from Bio import AlignIO, SeqIO, SeqUtils
term_size = get_terminal_size()
#AminoAcid = Literal["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "-"]
#Clustal = Literal['','.',':','*']
#SeqClass=new_class(SeqRecord)
#Position = Tuple[int, Clustal, AminoAcid, AminoAcid, SeqClass, SeqClass]
#PairDiff = list[Position]
def format_print(input1: dict, input2: dict, header1: str, header2: str, modestring: str):
    print(modestring)
    print(header1)
    pprint(input1)
    print(header2)
    pprint(input2)
    print(u'\u2500' * term_size.columns)

def transform(p_list, pos: int) -> dict:
    test_a=Counter(list(list(zip(*p_list))[pos]))
    three_counter={}
    for key in test_a:
        three_counter[SeqUtils.seq3(key)] = test_a[key]
    if 'Xaa' in three_counter:
        three_counter['GAP'] = three_counter.pop('Xaa')
    return three_counter


parser = ArgumentParser(prog="get_diff",
description="a script that takes paired fasta sequences from two species and returns count of mismatched AAs"
                        )
subcommand = parser.add_mutually_exclusive_group(required=True)
subcommand.add_argument("--compare-combined", action="store_true", help="compare all proteins")
subcommand.add_argument("--compare-protein", action="store", type=int, help="index of protein to compare", metavar='i', nargs="?")
parser.add_argument("in_species", action="extend", nargs="+", help="species base used as input for getPairs.sh")
parser.add_argument("-a, --all", action="store_true", dest="all")
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
        pairlist = []
        for kndex, letter in enumerate(alignment.column_annotations['clustal_consensus']):
            name1=self.seqnamelist[i].description
            name2=self.seqidlist[i].description
            base1 = alignment[0, kndex]
            base2 = alignment[1, kndex]
            diff = letter
            loc = kndex
            if diff != "*":
                pairlist.append((loc, diff, base1, base2, name1, name2))
        return pairlist
if parsed.compare_protein:
    local = SpeciesCompare(parsed.in_species[0])
    REMAINING = 0
    parseList = []
    if parsed.all:
        REMAINING = len(local.seqidlist) - 1
        while REMAINING >1:
            parseList.append(local.getpairlist(REMAINING))
            REMAINING -= 1
    else:
        parseList.append(local.getpairlist(parsed.compare_protein-1))
    for parse in parseList:
        if parse != []:
            transA=transform(parse, 2)
            transB=transform(parse, 3)
            headA=parse[0][4]
            headB=parse[0][5]
            format_print(transA, transB, headA, headB, "")

        




if parsed.compare_combined:
    listComparisons: list[SpeciesCompare] = []
    for comparator in parsed.in_species:
        local = SpeciesCompare(comparator)
        testAll=[local.getpairlist(i-1) for i in range(len(local.seqidlist)+1)]
        flat=[item for sublist in testAll for item in sublist]
        dictA=transform(flat,2)
        dictB=transform(flat,3)
        filename=comparator.rsplit('/',1)[1]
        headerA = filename.split('v')[0] ##this is imprecise and probably a bug waiting to explode
        headerB = filename.split('v')[1]
        format_print(dictA, dictB, headerA, headerB, "Printing amino-acid differences for all proteins combined:")
