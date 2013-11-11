#!/bin/python
## Genetic stability tool, v. 0.1 Alpha Americano.
## Gets you DNA sequences which are most unlike the original, but coding for
## the same protein.
## Go to http://2011.igem.org/Team:Edinburgh/Genetic_instability .
## Copyright (C) 2011 Lukasz Kopec, iGEM 2011 Team Edinburgh project
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## This tool includes bits of Biopython code, released under Biopython
## license. The license is available online: 
## http://www.biopython.org/DIST/LICENSE and in Bio folder.
## All Biopython code is (C) Copyright 2000 by Jeffrey Chang.
## All rights reserved.
## 
##
##
## iGEM 2011 Team Edinburgh project
## Takes a coding DNA sequence and returns one of sequences which are
## most unlike it (using simple distance metric), but code for the same
## protein.
## 
## This one is actually using the genetic algorithm for solving DNA-related
## problems! Cool, eh?

import sys
sys.path.append('./Bio')
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Data import CodonTable
from math import log10
import random, string, sys, os.path, argparse
global fasta_id
fasta_id = "Unknown protein"
standard_table = CodonTable.unambiguous_dna_by_id[1]
standard_table.forward_table.update(dict([('TAA', '*'), ('TGA', '*'), ('TAG', '*')]))

code = {'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'K': ['AAA', 'AAG'],
'N': ['AAT', 'AAC'], 'M': ['ATG'],
'D': ['GAT', 'GAC'], 'F': ['TTT', 'TTC'],
'C': ['TGT', 'TGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
'Q': ['CAA', 'CAG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
'E': ['GAA', 'GAG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'W': ['TGG'],
'H': ['CAT', 'CAC'], 'Y': ['TAT', 'TAC'],
'I': ['ATT', 'ATC', 'ATA'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
'*': ['TAA', 'TGA', 'TAG']}

# Codon Usage in E. coli
# from HENAUT and DANCHIN:Analysis and Predictions from Escherichia coli sequences.
# Escherichia coli and Salmonella, Vol. 2, Ch. 114:2047-2066, 1996, Neidhardt FC ed., ASM press,
# Washington, D.C.

codon_usage = {'A': {'GCT': 0.2754, 'GCC': 0.1614, 'GCA': 0.2401, 'GCG': 0.3230},
'R': {'AGA': 0.0062, 'AGG': 0.0029, 'CGT': 0.6425, 'CGC': 0.3297, 'CGA': 0.0107, 'CGG': 0.0080},
'N': {'AAT': 0.1725, 'AAC': 0.8275},
'D': {'GAT': 0.4605, 'GAC': 0.5395},
'C': {'TGT': 0.3885, 'TGC': 0.6115},
'E': {'GAA': 0.7535, 'GAG': 0.2465},
'Q': {'CAA': 0.1865, 'CAG': 0.8135},
'G': {'GGT': 0.5084, 'GGC': 0.4283, 'GGA': 0.0197, 'GGG': 0.0436},
'H': {'CAT': 0.2977, 'CAC': 0.7023},
'I': {'ATT': 0.3349, 'ATC': 0.6594, 'ATA': 0.0057},
'L': {'TTA': 0.0344, 'TTG': 0.0547, 'CTT': 0.0556, 'CTC': 0.0803, 'CTA': 0.0083, 'CTG': 0.7667},
'K': {'AAA': 0.7855, 'AAG': 0.2145},
'M': {'ATG': 1.000},
'F': {'TTT': 0.2908, 'TTC': 0.7092},
'P': {'CCT': 0.1123, 'CCC': 0.0163, 'CCA': 0.1525, 'CCG': 0.7189},
'S': {'TCT': 0.3241, 'TCC': 0.2656, 'TCA': 0.0479, 'TCG': 0.0739, 'AGT': 0.0452, 'AGC': 0.2433},
'T': {'ACT': 0.2908, 'ACC': 0.5360, 'ACA': 0.0467, 'ACG': 0.1265},
'W': {'TGG': 1.000},
'Y': {'TAT': 0.3523, 'TAC': 0.6477},
'V': {'GTT': 0.3977, 'GTC': 0.1345, 'GTA': 0.1997, 'GTG': 0.2681},
'*': {'TAA': 1/3.0, 'TGA': 1/3.0, 'TAG': 1/3.0}}

def str_comp(str1, str2):
    """ Takes two strings of equal length
    Returns simple edit distance
    """
    assert len(str1)==len(str2), 'Unequal length strings'
    diff=0;
    for i in range(len(str1)):
        if str1[i]!=str2[i]:
            diff=diff+1
    return diff

# from http://snippets.dzone.com/posts/show/732
def w_choice(lst):
    """ Weighted choice from list lst
    """
    probMass = sum(a[1] for a in lst)
    n = random.uniform(0, probMass)
    for item, weight in lst:
        if n < weight:
            break
        n = n - weight
    return item

def give_sample(amino):
    """ Gives a codon for amino acid _amino_, taking into account codon bias
    """
    probs=[]
    # Build a list of codons which are above CUTOFF value (set later)
    for codon in code[amino]:
        if codon_usage[amino][codon] > CUTOFF:
            probs.append((codon, codon_usage[amino][codon]))
    return w_choice(probs)

def get_individual(protein):
    """ Gives a number of DNA sequences coding for _protein_, taking into account
    codon bias
    """
    sequences=[]
    for i in range(NO_SEQUENCES):
        sequences.append(string.join([give_sample(amino) for amino in protein], '').upper())
    sequences=string.join(sequences, '')
    return sequences

def do_mutation(sequence):
    """ Perform a random mutation of a DNA sequence, changin one codon to an
    analogous codon.
    """
    r=random.randint(0, len(sequence)/3-1)
    codon=sequence[3*r:3*(r+1)]
    amino=standard_table.forward_table[codon]
    newcodon=random.choice(code[amino])
    result = sequence[:3*r]+newcodon+sequence[3*(r+1):]
    return result

def do_crossover(a,b):
    """ Perform a crossing-over of _a_ and _b_ DNA sequences
    """
    r=random.randint(0, len(a)/3)
    res1=a[:3*r]+b[3*r:]
    res2=b[:3*r]+a[3*r:]
    return (res1, res2)

def get_codon_bias(seq):
    """ Gets sum of log probabilities of codons for the whole sequence _seq_
    """
    codonBias = 0
    for i in range(len(seq)/3):
        codon = seq[i*3:(i+1)*3]
        amino = standard_table.forward_table[codon]
        codonBias += log10(codon_usage[amino][codon])
    return codonBias

def split_into_chunks(list, n):
    """ Splits _list_ into chunks of size _n_
    """
    return [list[i:i+n] for i in range(0, len(list), n)]

def cmp_all(lst):
    """ Compares all items in a list with each other
    """
    result = []
    for i in range(len(lst)):
        for j in range(i+1, len(lst)):
            result.append(str_comp(lst[i], lst[j]))
    return result

def fitness(sequence, source):
    """ Returns a fitness of _sequence_ as the sum of distance from _source_
    and the overhead in codon bias from _source_
    """
    global codonBiasSource
    sequences = split_into_chunks(sequence, len(source))
    sequences.append(source)
    distance = sum(cmp_all(sequences))
    codonBiasAvg = sum([get_codon_bias(seq) for seq in sequences])/NO_SEQUENCES
    return distance + codonBiasAvg

def do_gen_alg(babyBoomers, source):
    """ The main body of the Genetic Algorithm
    """
    rand = random.random()
    generationX = babyBoomers #copy the previous generation, and then we will
    # edit it.
    if rand < P_CROSSOVER:
        for i in range(NO_INDS/2):
            crossover1 = random.choice(babyBoomers)
            crossover2 = random.choice(babyBoomers)
            #choose two individuals for crossing-over and make sure
            #they are different
            while crossover1==crossover2:
                crossover2 = random.choice(babyBoomers)
            (res1,res2) = do_crossover(crossover1, crossover2)
            generationX.append(res1)
            generationX.append(res2)
    rand = random.random()
    if rand < P_MUTATION:
        for i in range(NO_INDS/2):
            index1 = random.randint(0, len(generationX)-1)
            generationX[index1] = do_mutation(generationX[index1])
    generationX.sort(key = lambda x: fitness(x, source), reverse=True)
    #that should introduce some more genetic variation
    #while making sure best solutions get through
    result=generationX[:NO_INDS/2]
    selectFrom = [(ind, fitness(ind, source)) for ind in generationX[NO_INDS/2:]]
    for i in range(NO_INDS-NO_INDS/2):
        result.append(w_choice(selectFrom))
    return result

def sanity_check(source, results):
    """ Checks if all of results code for the same protein as source
    """
    return all([str(Seq(source).transcribe().translate()) == str(Seq(seq).transcribe().translate()) for seq in results])

def initialise_ga(protein, source):
    """ Initialises the Genetic Algorithm
    """
    global codonBiasSource
    codonBiasSource = get_codon_bias(source)
    AdamAndEve = [get_individual(protein) for i in range(NO_INDS)]
    for i in range(NO_GENS):
        AdamAndEve = do_gen_alg(AdamAndEve, source)
    return AdamAndEve

def get_largest_distance(source):
    """ Takes a source DNA fragment and returns one of sequences which are
        most unlike it (using simple distance metric), but code for the same
        protein.
    """
    assert len(source)%3 == 0, 'The tool only accepts genetic code in multiples of 3 bases, sequence given: ('+str(len(source))+' bp) '+source
    source=string.upper(source)
    if ('u' in source) and ('t' in source):
        print "WARNING: The input has both U and T"
    source=string.replace(source, 'U', 'T')
    assert all([ch in 'CGTA' for ch in source]), 'The tool only accepts DNA/RNA (please use C, G, T/U, A), sequence given: ('+str(len(source))+' bp) '+source

    protein = str(Seq(source).transcribe().translate())
    if MODE == 'g':
        #run genetic algorithm
        result = initialise_ga(protein, source)
        best = split_into_chunks(result[0], len(source))
        assert sanity_check(source, best)
        for i in range(NO_SEQUENCES):
            print SeqRecord(Seq(best[i]),
                   id=fasta_id+" copy: "+str(i),
                   description="Generated by Genetic Algorithm").format('fasta')
            print ""
        print ""
        best.append(source)
        for i in range(len(best)):
            for j in range(i+1, len(best)):
                print "Comp result "+str(i)+" and "+str(j)+" : "+str((str_comp(best[i], best[j])+0.0)/len(best[i]))
    elif MODE == 'b':
        #run best codon
        bestresult=[]
        for i in range(len(protein)):
            sourceCodon = source[i*3:i*3+3]
            current = code[protein[i]]
            current = sorted(current, key=lambda x: str_comp(sourceCodon,x), reverse=True)
            bestresult.append(current[0])
        bestresult = string.join(bestresult,'')
        print SeqRecord(Seq(bestresult),
                   id=fasta_id+" copy",
                   description="Generated by Best Codon method").format('fasta')
    elif MODE == 'r':
        #run random codon
        random.seed()
        randresult = string.join([random.choice(code[amino]) for amino in protein], '')
        print SeqRecord(Seq(randresult),
                   id=fasta_id+" copy",
                   description="Generated by Random Codon method").format('fasta')

def do_FASTA_parse(path):
    """ Parses a FASTA format file, and sets the global fasta_id to the file's ID
    """
    global fasta_id
    for seq_record in SeqIO.parse(path, "fasta"):
        fasta_id = seq_record.id
        get_largest_distance(str(seq_record.seq))

parser = argparse.ArgumentParser(description='Genetic stability tool.')
parser.add_argument('-f', '--file', help = 'Take a DNA sequence from specified file in FASTA format')
parser.add_argument('-d', '--dna', help = 'Use the specified DNA sequence')
parser.add_argument('-n', '--noseq', help = 'Number of sequences to generate', type=int, default=2)
parser.add_argument('-g', '--gen', help = 'Number of generations of Genetic Algorithm', type=int, default=50)
parser.add_argument('-i', '--ind', help = 'Number of individuals in Genetic Algorithm', type=int, default=100)
parser.add_argument('-c', '--cross', help = 'Probability of crossing-over (0<c<1)', type=float, default=0.5)
parser.add_argument('-m', '--mutate', help = 'Probability of mutation (0<m<1)', type=float, default=0.3)
parser.add_argument('-o', '--cutoff', help = 'Cutoff point of inefficient E.coli codons (0<o<0.1)', type=float, default=0.05)
parser.add_argument('-a', '--action', help = 'The mode of action - (b)est codon, (g)enetic algorthm, (r)andom codon', choices='bgr', default='g')
args=parser.parse_args()

MODE=args.action
NO_SEQUENCES = args.noseq
CUTOFF = args.cutoff
NO_GENS = args.gen
NO_INDS = args.ind
P_CROSSOVER = args.cross
P_MUTATION = args.mutate
try:
    if ((args.file) and (args.dna)):
        print "Error: Please specify either a filename or a DNA sequence."
    elif args.file:
        if os.path.exists(args.file):
            do_FASTA_parse(args.file)
        else:
            print "File not found!"
    elif args.dna:
        get_largest_distance(args.dna)
    elif ((not args.file) and (not args.dna)):
        print "Please input the DNA sequence: "
        source=raw_input()
        get_largest_distance(source)
except Exception as e:
    print e
