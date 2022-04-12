""" Takes fasta sequences from two ordered fasta files, in which each pair is a known alignment.
Developed against Python 3.9.7, BioPython 1.78, and NumPy 1.19.5
Tested with Muscle v3.8.1551.
"""
import plac
#from argparse import ArgumentParser
from collections import Counter
from io import StringIO
from os import get_terminal_size
from pprint import pprint

#from typing import Literal, Tuple
#from types import new_class
from Bio import AlignIO, SeqIO, SeqUtils
from skbio import Protein
from skbio_parasail import SubstitutionMatrix
from skbio_parasail import global_pairwise_align_protein as gpap

import calculate_cons_for_clustal_protein as cons


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
    print(u'\u2500' * term_size)

def transform(p_list, pos: int) -> dict:
    test_a=Counter(list(list(zip(*p_list))[pos]))
    three_counter={}
    for key in test_a:
        three_counter[SeqUtils.seq3(key)] = test_a[key]
    if 'Xaa' in three_counter:
        three_counter['GAP'] = three_counter.pop('Xaa')
    return three_counter


class SpeciesCompare():
    def __init__(self,base: str) -> None:
        self.seq_ids=SeqIO.parse("%sbestIDs.fsa"  % (base), "fasta")
        self.seqnames=SeqIO.parse("%sbestnames.fsa"  % (base), "fasta")
        self.seqidlist=list(self.seq_ids)
        self.seqnamelist=list(self.seqnames)
    def getpairlist(self, i: int):
        seqnameprot = Protein(str(self.seqnamelist[i].seq))
        seqidprot = Protein(str(self.seqidlist[i].seq))
        alignment = gpap(
            seqnameprot, seqidprot, substitution_matrix=SubstitutionMatrix.from_name("blosum90"))
        msa = alignment[0]
        handler = StringIO()
        _ = msa.write(handler, format='clustal')
        clustal = cons.calculate_cons_for_clustal_protein(_.getvalue())

        alignment = AlignIO.read(StringIO(clustal.getvalue()), "clustal")
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

def compared(input):
    local = SpeciesCompare(input)
    testAll=[local.getpairlist(i-1) for i in range(len(local.seqidlist)+1)]
    flat=[item for sublist in testAll for item in sublist]
    dictA=transform(flat,2)
    dictB=transform(flat,3)
    filename=input
    #filename=comparator.rsplit('/',1)[1]
    headerA = filename.split('v')[0] ##this is imprecise and probably a bug waiting to explode
    headerB = filename.split('v')[1]
    format_print(dictA, dictB, headerA, headerB, "Printing amino-acid differences for all proteins combined:")

def per(input, index, all):
        local = SpeciesCompare(input)
        REMAINING = 0
        parseList = []
        if all:
            REMAINING = len(local.seqidlist) - 1
            while REMAINING >1:
                parseList.append(local.getpairlist(REMAINING))
                REMAINING -= 1
        else:
            parseList.append(local.getpairlist(index-1))
        for parse in parseList:
            if parse != []:
                transA=transform(parse, 2)
                transB=transform(parse, 3)
                headA=parse[0][4]
                headB=parse[0][5]
                format_print(transA, transB, headA, headB, "")

commands = 'combinedCompare', 'perProtein'

@plac.annotations(input=('NamevID', 'positional'))
def combinedCompare(input):
    "Compare amino acids across all proteins in paired name/ID FASTAs"
    return(compared(input))

@plac.annotations(input=('NamevID', 'positional'),
                  index=('protein index', 'positional'),
                  all=('iterate through all', 'flag', 'a'))
def perProtein(input: str, index: int, all: bool):
    "Compare amino acids for specific protein"
    return(per(input, index, all))

def __missing__(name):
    return ('Unknown option: %r' % name)

def __exit__(etype, exc, tb):
    "Will be called automatically at the end of the interpreter loop"
    if etype in (None, GeneratorExit):
        print('ok')

main = __import__(__name__)
if __name__ == '__main__':  
    try:
        term_sized = get_terminal_size()
        term_size = term_sized.columns
    except: 
        term_size= 80

    plac.call(main)
