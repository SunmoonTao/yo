from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict
import sys
sys.path.append('/Users/ltao/PycharmProjects/yo/tao')
from myprimer import RandomDNA, rc

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
    # convert primers to mips
    Mly1_site="GAGTC"
    Mly1_F='TATGAGTGTGGAGTCGTTGC'
    Mly1_R='GCTTCCTGATGAGTCCGATG'
    OM6_spacer="NNNAGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCTNNN"
    raw=Mly1_F+rc(fwd)+OM6_spacer+rev+rc(Mly1_R)
    if (Mly1_site in raw[15:-15]) or (rc(Mly1_site) in raw[15:-15]):
        raw = '/5Phos/' + rc(fwd) + OM6_spacer + rev
        #return 'MlyI in amplicon, not suitable for MlyI mips, try earI adaptors'
    else:
        return raw

def primers2ssmips(fwd,rev):
    # convert primers to mips
    HTspacer="AGATCGGAAGAGCACACGTCTGAACTCTTTCCCTACACGACGCTCTTCCGATCT"
    raw='/5Phos/'+rc(fwd)+OM6_spacer+rev
    return raw

def yolist(start,stop,step=1):
    # adjust the stop side to include
    stop = stop + 1
    return list(range(start,stop,step))


def string_diversity(list_of_string):
    cycs = dict()

    for index in range(len(list_of_string[0])):
        cycle_score = 0
        x = []
        for s in list_of_string:
            x.append(s[index])
        pA = x.count('A') / float(len(x))
        pC = x.count('C') / float(len(x))
        pRed = pA + pC
        pT = x.count('T') / float(len(x))
        pG = x.count('G') / float(len(x))
        pGreen = pT + pG
        if 0.1 <= pA <= 0.3 and 0.1 <= pC <= 0.3 and 0.1 <= pT <= 0.3 and 0.1 <= pG <= 0.3:
            cycle_score = 100
        elif 0.25 <= pA + pC <= 0.75:
            cycle_score = 75
        elif pA + pC == 1 or pA + pC == 0:
            cycle_score = 0
        else:
            cycle_score = 50

        cycs[index] = cycle_score
    return cycs


def breaknumbers(n):
    if n <= max_e * 2:
        a = int(n / 2)
        b = n - a

        return (a, b)
    else:
        print('not possible')


def Fibonacci(n):
    if n < 0:
        print("Incorrect input")
        # First Fibonacci number is 0
    elif n == 0:
        return 0
    # Second Fibonacci number is 1
    elif n == 1:
        return 1
    else:
        return Fibonacci(n - 1) + Fibonacci(n - 2)

def stagger_generator(primer,number_of_staggers,score_threshold):
    count=1
    while count>0:
        if number_of_staggers < 4:
            print ('more staggers needed')
            break
        else:
            combo=list()
            combo.append(primer)
            prefixes=list()
            for prefix in [str(RandomDNA(i)) for i in yolist(1,number_of_staggers)]:
                combo.append(prefix+primer)
                prefixes.append(prefix)
            score=string_diversity(combo)
            if list(score.values()).count(0)==0 and sum(score.values())/len(primer)>score_threshold:
                print (sum(score.values()),sum(score.values())/len(primer))
                print (score.values())
    #             print (prefixes)
#                 for i in combo:
#                     print (i)
                count=0
                return prefixes
