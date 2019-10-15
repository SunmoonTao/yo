from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict


def rc(str_seq):
    return str(Seq(str_seq).reverse_complement())


def xls2dic(file_path,sheetname):
    idh4x = pd.read_excel(file_path, sheet_name="Sheet1")
    idh4x=idh4x.to_dict(orient="index")
    return idh4x


def replace_str_index(text,index=0,replacement=''):
    return text[:index]+replacement+text[index+1:] 

# input a string of dna seq, return a list of mutated seqs of this dna with one mutated locus
def srandom_mutate(dna):
    out=list() 
    mutations=['A',"T",'C','G']
    for i,s in enumerate(dna):
        for mut in mutations:
            out.append(dna[:i]+mut+dna[i+1:])
    return set(out)

# input a string of dna seq, return a list of mutated seqs of this dna with one mutated locus

def mrandom_mutate(dna,mutation_numbers=1):
    muts=defaultdict(set)
    for r in range(mutation_numbers+1):
        mutations=['A',"T",'C','G']
        if r==0:
            muts[r].add(dna)
        else:
            for dna in muts[r-1]:
                for i in range(len(dna)):
                    for mut in mutations:
                        muts[r].add(dna[:i]+mut+dna[i+1:])
    return muts[mutation_numbers]

# input a string of dna seq, return a list of mutated seqs of this dna with one mutated locus
def random_mutate(dna):
    out=list() 
    mutations=['A',"T",'C','G']
    for i,s in enumerate(dna):
        for mut in mutations:
            out.append(dna[:i]+mut+dna[i+1:])
    return set(out)

def mrandom_mutate(dna,mutation_numbers=1):
    muts=defaultdict(set)
    for r in range(mutation_numbers+1):
        mutations=['A',"T",'C','G']
        if r==0:
            muts[r].add(dna)
        else:
            for dna in muts[r-1]:
                for i in range(len(dna)):
                    for mut in mutations:
                        muts[r].add(dna[:i]+mut+dna[i+1:])
    return muts[mutation_numbers]

def primers2mips(fwd,rev):
    Mly1_site="GAGTC"
    Mly1_F='TATGAGTGTGGAGTCGTTGC'
    Mly1_R='GCTTCCTGATGAGTCCGATG'
    OM6_spacer="NNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNN"
    raw=Mly1_F+rc(fwd)+OM6_spacer+rev+rc(Mly1_R)
    if (Mly1_site in raw[15:-15]) or (rc(Mly1_site) in raw[15:-15]):
        output='MlyI in amplicon, not suitable for MlyI mips, try earI adaptors'
    else:
        output = raw
    return (output)
