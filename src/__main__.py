""" Takes fasta sequences from two ordered fasta files, in which each pair is a known alignment.
"""
from collections import Counter
from io import StringIO
from os import get_terminal_size
from pprint import pprint
from typing import cast, Iterable, Literal, NamedTuple

from Bio import AlignIO, SeqIO, SeqUtils
from Bio.Data import IUPACData

from matplotlib.pyplot import show, close, autoscale, ion
from Bio.Align import MultipleSeqAlignment

import plac

from skbio import Protein

from skbio_parasail import SubstitutionMatrix
from skbio_parasail import global_pairwise_align_protein as gpap

import calculate_cons_for_clustal_protein as cons
import render

clust = Literal[' ', '', '.', ':', '*']
AminoAcid = Literal['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 
                    'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y', '-', 'B', 'X', 'Z', 'J', 'U', 'O'] 


class BasePair(NamedTuple):
    pos: int
    clust: clust
    amino1: AminoAcid
    amino2: AminoAcid
    name1: str
    name2: str


def format_print(input1: dict, input2: dict, header1: str, header2: str, modestring: str) -> None:
    print(modestring)
    print(header1)
    pprint(input1)
    print(header2)
    pprint(input2)
    print('\u2500' * term_size)


def transform(p_list, pos: int) -> dict:
    test_a = Counter(list(list(zip(*p_list))[pos]))
    three_counter = {}
    for key in test_a:
        three_counter[SeqUtils.seq3(key)] = test_a[key]
    if 'Xaa' in three_counter:
        three_counter['GAP'] = three_counter.pop('Xaa')
    return three_counter


class SpeciesCompare():
    seq_ids: Iterable[SeqIO.SeqRecord]
    seqnames: Iterable[SeqIO.SeqRecord]
    seqidlist: list[SeqIO.SeqRecord]
    seqnamelist: list[SeqIO.SeqRecord]
    
    def __init__(self, base: str) -> None:
        self.seq_ids = SeqIO.parse("%sbestIDs.fsa" % (base), "fasta")
        self.seqnames = SeqIO.parse("%sbestnames.fsa" % (base), "fasta")
        self.seqidlist = list(self.seq_ids)
        self.seqnamelist = list(self.seqnames)

    def getpairlist(self, i: int) -> list[BasePair]:
        seqnameprot = Protein(str(self.seqnamelist[i].seq))
        seqidprot = Protein(str(self.seqidlist[i].seq))
        alignment = gpap(
            seqnameprot, seqidprot, substitution_matrix=SubstitutionMatrix.from_name("blosum90"))
        msa = alignment[0]
        handler = StringIO()
        _ = msa.write(handler, format='clustal')
        clustal = cons.calculate_cons_for_clustal_protein(_.getvalue())
        alignmenta = AlignIO.read(StringIO(clustal.getvalue()), "clustal")
        assert isinstance(alignmenta, MultipleSeqAlignment)
        alignment = cast(MultipleSeqAlignment, alignmenta)
        pairlist = []
        for kndex, letter in enumerate(alignment.column_annotations['clustal_consensus']):
            name1a = self.seqnamelist[i].description
            name2a = self.seqidlist[i].description
            base1a = alignment[0, kndex]
            base2a = alignment[1, kndex]
            diffa = letter
            loc = kndex
            try:
                assert diffa in [' ', '', '.', ':', '*'], diffa
                diff = cast(clust, diffa)
                assert base1a in IUPACData.extended_protein_letters + '-' + 'X', base1a
                base1 = cast(AminoAcid, base1a)
                assert base2a in IUPACData.extended_protein_letters + '-' + 'X', base2a
                base2 = cast(AminoAcid, base2a)
                name1 = cast(str, name1a)
                name2 = cast(str, name2a)
            except AssertionError as a:
                print(a)
                raise
            else:
                tup = cast(BasePair, (loc, diff, base1, base2, name1, name2))
                pairlist.append(tup)
        return pairlist


def compared(inputa, same):
    local = SpeciesCompare(inputa)
    testAll = [local.getpairlist(i-1) for i in range(len(local.seqidlist)+1)]
    if same:
        flat = [item for sublist in testAll for item in sublist if item[1] == '*']
    else:
        flat = [item for sublist in testAll for item in sublist if item[1] != '*']

    dictA = transform(flat, 2)
    dictB = transform(flat, 3)
    filename = inputa
    # filename=comparator.rsplit('/',1)[1]
    headerA = filename.split('v')[0]  # this is imprecise and probably a bug waiting to explode
    headerB = filename.split('v')[1]
    format_print(dictA, dictB, headerA, headerB, "Printing amino-acids for all proteins combined:")
    render.renderCount(dictA, dictB, headerA, headerB)


def per(inputa, index, all_, residuals, returnList=False):
    index = int(index)
    local = SpeciesCompare(inputa)
    REMAINING = 0
    parseList = []
    renderList = []
    if all_:
        REMAINING = len(local.seqidlist) - 1
        while REMAINING >= 1:
            parseList.append(local.getpairlist(REMAINING))
            REMAINING -= 1
    else:
        parseList.append(local.getpairlist(index-1))
    for parse in parseList:
        if parse != []:
            transA = transform(parse, 2)
            transB = transform(parse, 3)
            headA = parse[0][4]
            headB = parse[0][5]
            if residuals:
                if all_:

                    netDict = {
                        k: abs(v - transB[k] if k in transB.keys() else v) for k, v in transA.items()}
                    renderList.append((netDict, headA + '\n' + headB))

                else:
                    netDict = {k: (v - transB[k] if k in transB.keys() else v) for k, v in transA.items()}
                    render.renderSingle(netDict, headA + '\n' + headB)
                    show()
            else:
                format_print(transA, transB, headA, headB, "")
                render.renderCount(transA, transB, str(headA), str(headB))
    if returnList:
        return renderList
    if renderList != []:
        figIns = render.singleheatmapper(renderList)
        figIns.rendermulti()
        show()


commands = 'combinedCompare', 'perProtein', 'multiChart'


@plac.annotations(infile=('NamevID', 'positional'),
                  same=('match', 'flag', 'm'))
def combinedCompare(infile: str, same: bool):
    "Compare amino acids across all proteins in paired name/ID FASTAs"
    return(compared(infile, same))


@plac.annotations(infile=('NamevID', 'positional'),
                  index=('protein index', 'positional'),
                  all_=('iterate through all', 'flag', 'a'),
                  residuals=('print protein net difference', 'flag', 'r'))
def perProtein(infile: str, index: int, all_: bool, residuals: bool):
    "Compare amino acids for specific protein"
    return(per(infile, index, all_, residuals))

@plac.annotations(infiles=('multi NamevID', 'positional', None, str))
def multiChart(*infiles: str):
    '''like perProtein -r -a, across multiple files'''
    presublist = []
    namelist = []
    for nvi in infiles:
        presublist.append(per(nvi, 2, True, True, True))
        namelist.append(nvi)
    multimapper = render.multiheatmapper(presublist, namelist)
    for i in range(len(presublist)):
        multimapper.seltmp(i)
        multimapper.rendermulti(cbar=True if i + 1 == len(presublist) else False)
    close(2)
    autoscale(True, 'both', True)
    show(block=True)

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
    except ValueError:
        term_size = 80

    plac.call(main)
